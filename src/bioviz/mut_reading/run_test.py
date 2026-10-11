import os
from bioviz.mut_reading import MutationParser, StructuralCoordinateMapper


def test_mutation_pipeline():
    print("=== Testing Mutation Parsing & 3D Structural Mapping ===")

    # 1. Load mutations using default predetermined values or a custom list
    records = MutationParser.load_mutations()

    # 2. Check if the local structure file exists in the predetermined_proteins folder
    structure_path = "src/bioviz/mut_reading/predetermined_proteins/1crn.cif"
    if os.path.exists(structure_path):
        print(f"\nFound structure file: {structure_path}. Initializing StructuralCoordinateMapper...")

        # Initialize mapper (Crambin 1crn has chain A; passing a sample offset if needed)
        mapper = StructuralCoordinateMapper(structure_file=structure_path, chain_id="A", residue_offset=0)

        # Map mutation records to 3D spatial coordinates
        mapped_data = mapper.map_mutations(records)

        print("\n--- 3D Mapping Results Output ---")
        for res in mapped_data:
            print(
                f"Residue {res['structure_position']} ({res['wt_residue']}->{res['mut_residue']}) | Status: {res['status']}")
            if res['ca_coordinate']:
                print(f"  -> Alpha Carbon (CA) 3D Coords: {res['ca_coordinate']}")
            else:
                print(f"  -> Coordinates unavailable for this index.")
    else:
        print(f"\n[INFO] Structure file {structure_path} not found. Skipping 3D coordinate mapping test.")


if __name__ == "__main__":
    test_mutation_pipeline()