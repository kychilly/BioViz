import os
from typing import List, Dict, Any, Optional
from Bio import PDB
from .parser import MutationRecord


class StructuralCoordinateMapper:
    """Maps 1D mutation sequence positions to 3D structural atoms from PDB/CIF files."""

    def __init__(self, structure_file: str, chain_id: str = "A", residue_offset: int = 0):
        self.structure_file = structure_file
        self.chain_id = chain_id
        self.residue_offset = residue_offset
        self.structure = self._load_structure()

    def _load_structure(self) -> PDB.Structure.Structure:
        """Loads a PDB or CIF structure file using Biopython."""
        ext = os.path.splitext(self.structure_file)[1].lower()
        if ext == '.cif':
            parser = PDB.MMCIFParser(QUIET=True)
        elif ext == '.pdb':
            parser = PDB.PDBParser(QUIET=True)
        else:
            raise ValueError(f"Unsupported structure file format: {ext}")

        return parser.get_structure("protein_structure", self.structure_file)

    def map_mutations(self, mutations: List[MutationRecord]) -> List[Dict[str, Any]]:
        """
        Maps mutation positions to 3D coordinates, verifying wild-type residues,
        handling chain filtering, and applying residue index offsets.
        """
        model = self.structure[0]  # Default to the first model
        if self.chain_id not in model:
            raise ValueError(f"Chain '{self.chain_id}' not found in structure {self.structure_file}")

        chain = model[self.chain_id]
        mapped_results = []

        for mut in mutations:
            # Apply offset if the crystal structure numbering differs from standard sequence numbering
            target_resi = mut.position + self.residue_offset
            residue_id = (' ', target_resi, ' ')  # standard PDB residue ID tuple (hetfield, resseq, icode)

            mapping_status = "SUCCESS"
            ca_coord = None
            heavy_atoms = {}
            struct_wt = None

            if residue_id in chain:
                residue = chain[residue_id]
                struct_wt = residue.get_resname()

                # Convert 3-letter to 1-letter amino acid code check if possible
                from Bio.PDB.Polypeptide import three_to_one
                try:
                    one_letter_wt = three_to_one(struct_wt)
                    if one_letter_wt != mut.wt_residue:
                        mapping_status = f"WARNING_WT_MISMATCH(Struct:{one_letter_wt} vs MutFile:{mut.wt_residue})"
                except Exception:
                    pass

                # Extract Alpha Carbon (CA) coordinate
                if "CA" in residue:
                    ca_coord = residue["CA"].get_coord().tolist()

                # Extract sidechain heavy atom coordinates
                for atom in residue:
                    if atom.get_name() != "N" and atom.get_name() != "C" and atom.get_name() != "O":
                        heavy_atoms[atom.get_name()] = atom.get_coord().tolist()
            else:
                mapping_status = "ERROR_RESIDUE_NOT_FOUND"

            mapped_results.append({
                "protein_id": mut.protein_id,
                "original_position": mut.position,
                "structure_position": target_resi,
                "wt_residue": mut.wt_residue,
                "mut_residue": mut.mut_residue,
                "structural_wt": struct_wt,
                "chain": self.chain_id,
                "status": mapping_status,
                "ca_coordinate": ca_coord,
                "sidechain_atoms": heavy_atoms
            })

        return mapped_results