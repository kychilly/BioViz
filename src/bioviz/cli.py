import argparse
import os
from .parser import MoleculeParser
from .engine import PyMoLEngine
from .compiler import VideoCompiler


def main():
    parser = argparse.ArgumentParser(description="BioViz Advanced Multi-State Binding Simulation Pipeline.")
    parser.add_argument("--apo_id", type=str, default="4ake",
                        help="PDB ID for open/unbound state (default: Adenylate Kinase 4AKE).")
    parser.add_argument("--holo_id", type=str, default="1ake",
                        help="PDB ID for closed/bound state (default: Adenylate Kinase 1AKE).")
    parser.add_argument("--output", type=str, default="binding_simulation.mp4", help="Path for the final output video.")
    parser.add_argument("--frames", type=int, default=60, help="Number of frames for the binding sequence.")

    args = parser.parse_args()

    print(f"=== Step 1: Acquiring & Aligning Multi-State Data ({args.apo_id.upper()} vs {args.holo_id.upper()}) ===")
    paths = MoleculeParser.fetch_multi_state_complex(apo_id=args.apo_id, holo_id=args.holo_id)

    apo_protein = paths["apo_protein"]
    holo_protein = paths["holo_protein"]
    ligand_path = paths["ligand"]

    # Align structures to anchor backbones
    aligned_apo = os.path.join("sample_data", "protein_apo_aligned.pdb")
    MoleculeParser.align_structures(mobile_pdb=apo_protein, target_pdb=holo_protein, output_aligned_path=aligned_apo)

    print("=== Step 3-7: Initializing Engine & Rendering Binding Trajectory ===")
    output_frames_dir = "rendered_frames"
    engine = PyMoLEngine(headless=True)
    engine.render_full_binding_simulation(
        apo_protein=aligned_apo,
        holo_protein=holo_protein,
        ligand_path=ligand_path,
        output_frames_dir=output_frames_dir,
        num_frames=args.frames
    )

    print("=== Step 8: Compiling Frames into Final MP4 Video ===")
    VideoCompiler.create_video_from_frames(output_frames_dir, args.output, fps=20)
    print(f"\n🎉 Multi-state binding simulation successfully compiled at: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()