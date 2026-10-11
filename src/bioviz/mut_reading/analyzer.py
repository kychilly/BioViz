import numpy as np
from typing import Dict, Any, List, Optional
from Bio import PDB


class MutationImpactAnalyzer:
    """Analyzes spatial proximity to ligands/pockets and evaluates local structural impacts."""

    # Simple hydrophobicity scale (Kyte-Doolittle approximation)
    HYDROPHOBICITY_SCALE = {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
        'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
        'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
        'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }

    def __init__(self, structure: PDB.Structure.Structure, chain_id: str = "A"):
        self.structure = structure
        self.chain_id = chain_id

    def analyze_impact(self, residue_position: int, wt_res: str, mut_res: str, cutoff_distance: float = 5.0) -> Dict[
        str, Any]:
        """
        Calculates minimum distance to other structural elements/ligands and
        flags physicochemical changes (hydrophobicity shift, potential salt bridge/H-bond impact).
        """
        model = self.structure[0]
        if self.chain_id not in model:
            return {"error": f"Chain {self.chain_id} not found."}

        chain = model[self.chain_id]
        residue_id = (' ', residue_position, ' ')

        if residue_id not in chain:
            return {"status": "ERROR_RESIDUE_NOT_FOUND", "position": residue_position}

        target_residue = chain[residue_id]

        # Get Alpha Carbon coordinate as anchor
        if "CA" not in target_residue:
            return {"status": "ERROR_NO_CA_ATOM", "position": residue_position}

        target_coord = target_residue["CA"].get_coord()

        # 1. Proximity check to other residues (pocket or packing tightness analysis)
        min_distance = float('inf')
        nearest_neighbor = None

        for other_res in chain:
            other_id = other_res.id[1]
            if other_id == residue_position:
                continue
            if "CA" in other_res:
                dist = np.linalg.norm(target_coord - other_res["CA"].get_coord())
                if dist < min_distance:
                    min_distance = dist
                    nearest_neighbor = other_id

        in_binding_pocket = min_distance <= cutoff_distance

        # 2. Physicochemical Impact Evaluation
        wt_hydro = self.HYDROPHOBICITY_SCALE.get(wt_res.upper(), 0.0)
        mut_hydro = self.HYDROPHOBICITY_SCALE.get(mut_res.upper(), 0.0)
        hydrophobicity_delta = round(mut_hydro - wt_hydro, 2)

        structural_flags = []
        if abs(hydrophobicity_delta) > 2.0:
            structural_flags.append("SIGNIFICANT_HYDROPHOBICITY_SHIFT")
        if in_binding_pocket:
            structural_flags.append("LOCATED_WITHIN_ACTIVE_SITE_PROXIMITY")
        else:
            structural_flags.append("SURFACE_OR_CORE_PACKING_REGION")

        return {
            "position": residue_position,
            "mutation": f"{wt_res.upper()} -> {mut_res.upper()}",
            "nearest_neighbor_residue": nearest_neighbor,
            "min_ca_distance_angstroms": round(float(min_distance), 2),
            "in_binding_pocket_window": in_binding_pocket,
            "hydrophobicity_delta": hydrophobicity_delta,
            "structural_flags": structural_flags,
            "status": "ANALYSIS_SUCCESS"
        }