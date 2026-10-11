import os
import json
import csv
import numpy as np
from typing import Dict, Any, List, Tuple
from Bio import PDB


class MutationReporter:
    """Handles serialization of modified structures and compilation of JSON/CSV summary reports."""

    @staticmethod
    def _convert_numpy_types(obj: Any) -> Any:
        """Recursively converts NumPy data types to native Python types for JSON serialization."""
        if isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, dict):
            return {k: MutationReporter._convert_numpy_types(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [MutationReporter._convert_numpy_types(item) for item in obj]
        return obj

    @staticmethod
    def save_structure(structure: PDB.Structure.Structure, output_pdb_path: str) -> str:
        """Saves a BioPython structure object to a specified PDB file path."""
        out_dir = os.path.dirname(output_pdb_path)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(output_pdb_path)
        print(f"[SUCCESS] Serialized modified structure to: {output_pdb_path}")
        return output_pdb_path

    @staticmethod
    def compile_reports(
            mutation_records: List[Dict[str, Any]],
            output_dir: str = "output",
            base_filename: str = "mutation_analysis_report"
    ) -> Tuple[str, str]:
        """
        Compiles comprehensive summary reports into JSON and CSV formats for downstream pipelines.
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

        json_path = os.path.join(output_dir, f"{base_filename}.json")
        csv_path = os.path.join(output_dir, f"{base_filename}.csv")

        # Sanitize and convert records to handle NumPy boolean/float types
        sanitized_records = [MutationReporter._convert_numpy_types(record) for record in mutation_records]

        # 1. Save JSON Report
        report_data = {
            "summary_metadata": {
                "total_mutations_analyzed": len(sanitized_records),
                "pipeline_version": "1.0.0"
            },
            "mutations": sanitized_records
        }
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(report_data, f, indent=4)

        # 2. Save CSV Report
        csv_columns = [
            "position", "mutation", "status", "clash_score",
            "min_ca_distance_angstroms", "hydrophobicity_delta", "structural_flags"
        ]
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=csv_columns, extrasaction='ignore')
            writer.writeheader()
            for record in sanitized_records:
                row = record.copy()
                if isinstance(row.get("structural_flags"), list):
                    row["structural_flags"] = ", ".join(row["structural_flags"])
                writer.writerow(row)

        print(f"[SUCCESS] Reports compiled successfully:\n  -> JSON: {json_path}\n  -> CSV: {csv_path}")
        return json_path, csv_path