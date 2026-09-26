import os
import pymol


class PyMoLEngine:
    def __init__(self, headless: bool = True):
        # Initialize PyMOL (-qc means quiet mode and no GUI)
        args = ['pymol', '-qc'] if headless else ['pymol']
        pymol.finish_launching(args)
        self.session_initialized = True

    def setup_scene(self, protein_path: str, ligand_path: str, output_dir: str):
        """Loads structures, styles them, and prepares for rendering frames."""
        os.makedirs(output_dir, exist_ok=True)

        # Clear any existing state
        pymol.cmd.reinitialize()

        # Load files
        pymol.cmd.load(protein_path, "protein")
        pymol.cmd.load(ligand_path, "ligand")

        # Set visual representations
        pymol.cmd.hide("everything", "all")
        pymol.cmd.show("cartoon", "protein")
        pymol.cmd.color("slate", "protein")

        pymol.cmd.show("sticks", "ligand")
        pymol.cmd.color("cyan", "ligand")

        # Orient view around the ligand active site
        pymol.cmd.zoom("ligand", buffer=5.0)
        pymol.cmd.bg_color("white")

    def render_trajectory_frames(self, output_dir: str, num_frames: int = 60):
        """Renders sequential frames (e.g., simulating a docking rotation or state transition)."""
        for i in range(num_frames):
            # Apply a slight rotation to simulate binding motion or camera sweep
            pymol.cmd.turn("y", 360.0 / num_frames)

            frame_path = os.path.join(output_dir, f"frame_{i:04d}.png")
            # Ray-trace for high quality, then save PNG
            pymol.cmd.ray(1200, 800)
            pymol.cmd.png(frame_path)

        return output_dir