import argparse
import os
import urllib.request
import pymol
from .parser import MoleculeParser
from .engine import PyMoLEngine
from .compiler import VideoCompiler


def main():
    parser = argparse.ArgumentParser(
        description="Automated PyMOL video rendering pipeline for protein-ligand binding simulations.")
    parser.add_argument("--pdb_id", type=str, default="1iep", help="RCSB PDB ID to download and render.")
    parser.add_argument("--output", type=str, default="binding_simulation.mp4", help="Path for the final output video.")
    parser.add_argument("--frames", type=int, default=45, help="Number of frames in the rotation loop.")

    args = parser.parse_args()

    print(f"=== Step 1: Fetching PDB: {args.pdb_id.upper()} ===")
    data_dir = "sample_data"
    os.makedirs(data_dir, exist_ok=True)

    raw_pdb_path = os.path.join(data_dir, f"{args.pdb_id}.pdb")
    protein_path = os.path.join(data_dir, "protein.pdb")
    ligand_path = os.path.join(data_dir, "ligand.pdb")

    if not os.path.exists(raw_pdb_path):
        url = f"https://files.rcsb.org/download/{args.pdb_id}.pdb"
        print(f"Downloading from {url}...")
        urllib.request.urlretrieve(url, raw_pdb_path)

    validated_path = MoleculeParser.validate_file_path(raw_pdb_path)

    print("=== Step 2: Preparing Coordinates ===")
    pymol.finish_launching(['pymol', '-qc'])
    pymol.cmd.load(validated_path, "complex")
    pymol.cmd.select("prot_sel", "polymer")
    pymol.cmd.select("lig_sel", "organic")
    pymol.cmd.save(protein_path, "prot_sel")
    pymol.cmd.save(ligand_path, "lig_sel")
    pymol.cmd.delete("all")

    print("=== Step 3: Rendering Frames ===")
    output_frames_dir = "rendered_frames"
    engine = PyMoLEngine(headless=True)
    engine.setup_scene(protein_path, ligand_path, output_frames_dir)
    engine.render_trajectory_frames(output_frames_dir, num_frames=args.frames)

    print("=== Step 4: Compiling Video ===")
    VideoCompiler.create_video_from_frames(output_frames_dir, args.output, fps=15)
    print(f"\n🎉 Simulation video successfully generated at: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()