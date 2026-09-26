import os

class MoleculeParser:
    @staticmethod
    def validate_file_path(file_path: str) -> str:
        """Validates that a given molecular coordinate file exists."""
        abs_path = os.path.abspath(file_path)
        if not os.path.exists(abs_path):
            raise FileNotFoundError(f"Molecular file not found at: {abs_path}")
        return abs_path