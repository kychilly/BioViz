import os
import pymol
import numpy as np
from bioviz.rendering.parser import MoleculeParser
from bioviz.rendering.compiler import VideoCompiler


class SimpleLigandApproachEngine:
    def __init__(self, headless: bool = True):
        self.headless = headless
        try:
            args = ['pymol', '-qc'] if self.headless else ['pymol']
            pymol.finish_launching(args)
        except Exception:
            pass

    def render_ligand_approach(self, protein_path: str, ligand_path: str,
                               output_frames_dir: str, num_frames: int = 60):
        os.makedirs(output_frames_dir, exist_ok=True)
        pymol.cmd.delete("all")

        # Load protein and ligand
        pymol.cmd.load(protein_path, "protein")
        pymol.cmd.load(ligand_path, "ligand_template")

        # Style the protein
        pymol.cmd.show("cartoon", "protein")
        pymol.cmd.color("slate", "protein")
        pymol.cmd.bg_color("white")
        pymol.cmd.set("ray_shadows", 0)
        pymol.cmd.set("antialias", 2)

        base_coords = []
        pymol.cmd.iterate_state(1, "ligand_template", "base_coords.append((x, y, z))",
                                space={'base_coords': base_coords})
        base_coords = np.array(base_coords, dtype=float)

        pymol.cmd.create("ligand", "ligand_template", 1, 1)
        for f in range(2, num_frames + 1):
            pymol.cmd.create("ligand", "ligand_template", 1, f)
        pymol.cmd.delete("ligand_template")

        # Style the ligand
        pymol.cmd.show("sticks", "ligand")
        pymol.cmd.show("spheres", "ligand")
        pymol.cmd.set("sphere_scale", 0.25, "ligand")
        pymol.cmd.color("cyan", "ligand")

        # Focus camera on the active site / ligand region
        pymol.cmd.frame(1)
        pymol.cmd.zoom("protein within 6 of ligand", buffer=3.0)
        locked_view = pymol.cmd.get_view()

        entry_offset = np.array([0.0, 0.0, 25.0])

        print(f"Rendering ligand approach animation across {num_frames} frames...")

        for i in range(1, num_frames + 1):
            fraction = (i - 1) / max(1, num_frames - 1)
            pymol.cmd.frame(i)

            # Fraction goes from 1.0 (far away) down to 0.0 (fully docked in active site)
            t_factor = 1.0 - fraction
            current_shift = entry_offset * t_factor
            shifted_coords = base_coords + current_shift

            l_idx = [0]

            def update_ligand_coords(x, y, z):
                idx = l_idx[0]
                if idx < len(shifted_coords):
                    val = shifted_coords[idx]
                    l_idx[0] += 1
                    return list(val)
                return [x, y, z]

            # Update ligand state for the current frame
            pymol.cmd.alter_state(
                i,
                "ligand",
                "(x, y, z) = update_ligand_coords(x, y, z)",
                space={'update_ligand_coords': update_ligand_coords}
            )

            # Lock the camera view so the scene doesn't move
            pymol.cmd.set_view(locked_view)

            # Save frame
            frame_path = os.path.join(output_frames_dir, f"frame_{i:04d}.png")
            pymol.cmd.png(frame_path, width=1200, height=900, ray=0)
            print(f"Rendered frame {i}/{num_frames}")

        pymol.cmd.delete("all")
        print("Frame rendering complete.")
        return output_frames_dir


def main():
    pdb_id = "1ake"
    print(f"Fetching structure {pdb_id.upper()}...")
    paths = MoleculeParser.process_arbitrary_pdb(pdb_id=pdb_id)

    protein_file = paths["protein"]
    ligand_file = paths["ligand"]

    if not ligand_file:
        raise ValueError(f"No ligand found in {pdb_id.upper()}. Please choose a PDB ID with a bound ligand.")

    output_frames_dir = "ligand_approach_frames"
    engine = SimpleLigandApproachEngine(headless=True)
    engine.render_ligand_approach(
        protein_path=protein_file,
        ligand_path=ligand_file,
        output_frames_dir=output_frames_dir,
        num_frames=20
    )

    output_dir = "rendering_videos"
    os.makedirs(output_dir, exist_ok=True)
    output_video = os.path.join(output_dir, "ligand_approach.mp4")

    print("Compiling video...")
    VideoCompiler.create_video_from_frames(output_frames_dir, output_video, fps=5)
    print(f"\n🎉 Video successfully generated at: {os.path.abspath(output_video)}")


if __name__ == "__main__":
    main()