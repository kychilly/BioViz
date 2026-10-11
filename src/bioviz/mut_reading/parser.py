import csv
import json
import os
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, field


@dataclass
class MutationRecord:
    """Standardized data container for a single mutation entry."""
    protein_id: str
    position: int
    wt_residue: str
    mut_residue: str
    chrom: Optional[str] = None
    additional_info: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        self.wt_residue = self.wt_residue.upper()
        self.mut_residue = self.mut_residue.upper()


class MutationParser:
    """Parser utility to ingest mutation data from standard bioinformatics formats."""

    @staticmethod
    def parse_csv_tsv(filepath: str, delimiter: str = ',') -> List[MutationRecord]:
        """Parses comma-separated or tab-separated tables containing mutations."""
        records = []
        with open(filepath, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            for row in reader:
                protein_id = row.get("PROTEIN_ID") or row.get("CHROM") or "UNKNOWN"
                pos = int(row.get("POS") or row.get("POSITION") or 0)
                wt = row.get("WT") or row.get("REF") or ""
                mut = row.get("MUT") or row.get("ALT") or ""

                records.append(MutationRecord(
                    protein_id=protein_id,
                    position=pos,
                    wt_residue=wt,
                    mut_residue=mut,
                    chrom=row.get("CHROM"),
                    additional_info=row
                ))
        return records

    @staticmethod
    def parse_json(filepath: str) -> List[MutationRecord]:
        """Parses structured JSON mutation profiles."""
        records = []
        with open(filepath, mode='r', encoding='utf-8') as f:
            data = json.load(f)
            mutations_list = data if isinstance(data, list) else data.get("mutations", [])
            for item in mutations_list:
                records.append(MutationRecord(
                    protein_id=item.get("protein_id") or item.get("chrom", "UNKNOWN"),
                    position=int(item.get("position", 0)),
                    wt_residue=str(item.get("wt_residue", item.get("ref", ""))),
                    mut_residue=str(item.get("mut_residue", item.get("alt", ""))),
                    chrom=item.get("chrom"),
                    additional_info=item
                ))
        return records

    @classmethod
    def load_mutations(cls, filepath: Optional[str] = None, file_type: Optional[str] = None) -> List[MutationRecord]:
        """
        Loads mutations from a specified file and format.
        If no filepath is provided, it falls back to a default predetermined built-in source.
        """
        if not filepath:
            print("[INFO] No mutation file path provided. Falling back to default predetermined dataset...")
            return cls._get_default_predetermined_mutations()

        # If file_type is explicitly provided, use it; otherwise infer from extension
        ext = f".{file_type.lower().strip('.')}" if file_type else os.path.splitext(filepath)[1].lower()

        if ext == '.csv':
            return cls.parse_csv_tsv(filepath, delimiter=',')
        elif ext in ['.tsv', '.txt']:
            return cls.parse_csv_tsv(filepath, delimiter='\t')
        elif ext == '.json':
            return cls.parse_json(filepath)
        else:
            raise ValueError(f"Unsupported or unrecognized file type/extension: {ext}")

    @staticmethod
    def _get_default_predetermined_mutations() -> List[MutationRecord]:
        """Fallback predetermined mutations if no file path is supplied."""
        default_data = [
            {"protein_id": "EGFR", "position": 858, "wt_residue": "L", "mut_residue": "R"},
            {"protein_id": "TEM1", "position": 104, "wt_residue": "E", "mut_residue": "K"},
            {"protein_id": "EGFR", "position": 790, "wt_residue": "T", "mut_residue": "M"}
        ]
        return [
            MutationRecord(
                protein_id=item["protein_id"],
                position=item["position"],
                wt_residue=item["wt_residue"],
                mut_residue=item["mut_residue"]
            ) for item in default_data
        ]