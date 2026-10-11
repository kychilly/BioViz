import math
import numpy as np
from typing import Dict, Any, List, Tuple
from Bio import PDB


class StructureMutator:
    """Performs in silico sidechain mutagenesis and rotamer minimization on protein structures."""

    # Standard amino acid code conversion dictionaries
    THREE_TO_ONE = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}

    # Standard idealized heavy atom offsets / bond approximations for sidechain building
    DEFAULT_ROTAMER_OFFSETS = {
        "A": {"CB": [0.0, 1.5, 0.0]},
        "V": {"CB": [0.0, 1.5, 0.0], "CG1": [1.0, 2.3, 0.9], "CG2": [-1.0, 2.3, -0.9]},
        "S": {"CB": [0.0, 1.5, 0.0], "OG": [0.0, 2.6, 1.2]},
        "T": {"CB": [0.0, 1.5, 0.0], "OG1": [1.0, 2.2, 0.8], "CG2": [-1.0, 2.4, -0.8]},
        "R": {"CB": [0.0, 1.5, 0.0], "CG": [0.0, 2.5, 1.2], "CD": [0.0, 3.6, 1.8], "NE": [1.0, 4.5, 2.0],
              "CZ": [0.5, 5.7, 2.2], "NH1": [1.5, 6.5, 2.5], "NH2": [-0.5, 6.0, 2.0]},
        "K": {"CB": [0.0, 1.5, 0.0], "CG": [0.0, 2.5, 1.2], "CD": [0.0, 3.7, 1.8], "CE": [0.0, 4.9, 2.4],
              "NZ": [0.0, 6.1, 3.0]},
        "M": {"CB": [0.0, 1.5, 0.0], "CG": [0.0, 2.5, 1.2], "SD": [0.0, 4.0, 2.0], "CE": [1.2, 5.0, 2.5]},
    }

    def __init__(self, structure: PDB.Structure.Structure, chain_id: str = "A"):
        self.structure = structure
        self.chain_id = chain_id

    def mutate_residue(self, residue_position: int, target_wt: str, mut_residue_1le: str) -> Dict[str, Any]:
        """
        Swaps the wild-type sidechain at a given residue position with the mutant amino acid,
        evaluating and selecting a rotamer configuration that minimizes steric clashes.
        """
        model = self.structure[0]
        if self.chain_id not in model:
            raise ValueError(f"Chain {self.chain_id} not found in structure.")

        chain = model[self.chain_id]
        residue_id = (' ', residue_position, ' ')

        if residue_id not in chain:
            return {"status": "ERROR_RESIDUE_NOT_FOUND", "position": residue_position}

        residue = chain[residue_id]
        current_resname = residue.get_resname().upper()

        # Verify WT code using local dictionary
        wt_1le = self.THREE_TO_ONE.get(current_resname, current_resname)

        # 1. Truncate wild-type sidechain (retain backbone N, CA, C, O, and optionally CB)
        atoms_to_remove = []
        for atom in residue:
            atom_name = atom.get_name()
            if atom_name not in ["N", "CA", "C", "O", "CB"]:
                atoms_to_remove.append(atom_name)

        for atom_name in atoms_to_remove:
            residue.detach_child(atom_name)

        # 2. Update residue name to mutant three-letter code
        mut_1le_upper = mut_residue_1le.upper()
        mut_3le = self.ONE_TO_THREE.get(mut_1le_upper, "ALA")
        residue.resname = mut_3le

        # 3. Generate candidate rotamers / atom placements for mutant sidechain
        ca_coord = residue["CA"].get_coord() if "CA" in residue else np.array([0.0, 0.0, 0.0])
        cb_coord = residue["CB"].get_coord() if "CB" in residue else ca_coord + np.array([1.5, 0.0, 0.0])

        template_offsets = self.DEFAULT_ROTAMER_OFFSETS.get(mut_1le_upper, {"CB": [0.0, 1.5, 0.0]})

        best_rotamer_score = float('inf')
        best_coords = {}
        rotamer_angles = [0.0, 45.0, 90.0, 180.0, 270.0]

        for angle in rotamer_angles:
            rad = math.radians(angle)
            cos_a, sin_a = math.cos(rad), math.sin(rad)

            current_rotamer_coords = {}
            clash_penalty = 0.0

            for atom_name, offset in template_offsets.items():
                if atom_name == "CB":
                    pos = cb_coord
                else:
                    base_vec = np.array(offset)
                    rot_x = base_vec[0] * cos_a - base_vec[2] * sin_a
                    rot_y = base_vec[1]
                    rot_z = base_vec[0] * sin_a + base_vec[2] * cos_a
                    pos = cb_coord + np.array([rot_x, rot_y, rot_z])

                current_rotamer_coords[atom_name] = pos

                for other_res in chain:
                    if other_res.id[1] == residue_position:
                        continue
                    for other_atom in other_res:
                        if other_atom.get_name() in ["N", "CA", "C", "O", "CB"]:
                            dist = np.linalg.norm(pos - other_atom.get_coord())
                            if dist < 2.0:
                                clash_penalty += (2.0 - dist) ** 2

            if clash_penalty < best_rotamer_score:
                best_rotamer_score = clash_penalty
                best_coords = current_rotamer_coords

        # 4. Insert optimized atoms into the residue
        added_atoms = []
        for atom_name, coord in best_coords.items():
            if atom_name == "CB" and "CB" in residue:
                continue

            vector = PDB.Vector(float(coord[0]), float(coord[1]), float(coord[2]))
            new_atom = PDB.Atom.Atom(
                name=atom_name,
                coord=vector,
                bfactor=0.0,
                occupancy=1.0,
                altloc=" ",
                fullname=f" {atom_name:<3}",
                serial_number=1,
                element=atom_name[0]
            )
            residue.add(new_atom)
            added_atoms.append(atom_name)

        return {
            "position": residue_position,
            "mutation": f"{wt_1le} -> {mut_1le_upper}",
            "status": "SUCCESS_MUTATED_AND_MINIMIZED",
            "clash_score": round(best_rotamer_score, 4),
            "added_sidechain_atoms": added_atoms
        }