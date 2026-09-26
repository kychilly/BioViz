import os
import urllib.request
import pymol
from bioviz import MoleculeParser, PyMoLEngine, VideoCompiler


def main():
    print("=== Step 1: Downloading Sample PDB (1IEP: Kinase-Inhibitor Complex) ===")
    pdb_id = "1iep"
    data_dir = "sample_data"
    os.makedirs(data_dir, exist_ok=True)

    raw_pdb_path = os.path.join(data_dir, f"{pdb_id}.pdb")
    protein_path = os.path.join(data_dir, "protein.pdb")
    ligand_path = os.path.join(data_dir, "ligand.pdb")

    # Download raw PDB if not already present
    if not os.path.exists(raw_pdb_path):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"Fetching {url}...")
        urllib.request.urlretrieve(url, raw_pdb_path)

    # Validate the file using your MoleculeParser package
    validated_path = MoleculeParser.validate_file_path(raw_pdb_path)

    # Use PyMOL headless to split protein and ligand into separate clean files
    print("=== Step 2: Preparing Protein and Ligand Coordinates ===")
    pymol.finish_launching(['pymol', '-qc'])
    pymol.cmd.load(validated_path, "complex")

    # Extract polymer (protein) and organic molecule (ligand)
    pymol.cmd.select("prot_sel", "polymer")
    pymol.cmd.select("lig_sel", "organic")

    pymol.cmd.save(protein_path, "prot_sel")
    pymol.cmd.save(ligand_path, "lig_sel")
    pymol.cmd.delete("all")

    print(alignment_msg := "Saved cleaned coordinates successfully.")

    print("=== Step 3: Initializing PyMoLEngine & Rendering Frames ===")
    output_frames_dir = "rendered_frames"
    output_video_path = "binding_simulation.mp4"

    # Initialize engine and setup scene
    engine = PyMoLEngine(headless=True)
    engine.setup_scene(protein_path, ligand_path, output_frames_dir)

    # Render 45 frames (smooth 360-degree rotation around the active site)
    engine.render_trajectory_frames(output_frames_dir, num_frames=45)
    print(f"Frames rendered to: {output_frames_dir}/")

    print("=== Step 4: Compiling Frames into MP4 Video ===")
    VideoCompiler.create_video_from_frames(output_frames_dir, output_video_path, fps=15)

    print(f"\n🎉 Test complete! Your video is ready at: {os.path.abspath(output_video_path)}")


if __name__ == "__main__":
    main()