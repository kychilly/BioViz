import os
import cv2


class StructuralIntegrator:
    """Translates variant residue positions into PyMOL scripts, renders local 3D structural images, and compiles 360 rotation videos."""

    def __init__(self, protein_object_name: str = "protein", output_dir: str = "src/bioviz/mutations/PyMOL_rendering"):
        self.protein_object_name = protein_object_name
        self.output_dir = output_dir

        # Ensure the rendering output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

    def generate_pymol_script(self, annotated_variants: list) -> str:
        """Generates PyMOL commands to highlight mutated residues and save an image."""
        affected_residues = set()
        for v in annotated_variants:
            pos = v["position"]
            aa_pos = (pos // 3) + 1
            affected_residues.add(str(aa_pos))

        if not affected_residues:
            return "# No structural mutations detected to highlight."

        resi_string = "+".join(sorted(affected_residues, key=int))

        # Use forward slashes for PyMOL path compatibility on Windows
        output_image_path = os.path.abspath(os.path.join(self.output_dir, "mutation_render.png")).replace("\\", "/")

        script_lines = [
            f"# Auto-generated PyMOL script for BioViz Mutation Analysis",
            f"cmd.hide('everything', '{self.protein_object_name}')",
            f"cmd.show('cartoon', '{self.protein_object_name}')",
            f"cmd.color('grey80', '{self.protein_object_name}')",
            f"",
            f"# Select and highlight mutated residues",
            f"cmd.select('mutated_sites', 'resi {resi_string} and {self.protein_object_name}')",
            f"cmd.show('sticks', 'mutated_sites')",
            f"cmd.color('firebrick', 'mutated_sites')",
            f"",
            f"# Set viewport dimensions, focus camera framing, and render image",
            f"cmd.viewport(1000, 1000)",
            f"cmd.zoom('mutated_sites', buffer=5.0)",
            f"cmd.ray(1000, 1000)",
            f"cmd.png('{output_image_path}')"
        ]

        return "\n".join(script_lines)

    def render_structure_image(self, annotated_variants: list, PDB_file_path: str = None) -> str:
        """
        Executes the rendering pipeline locally. If a PDB file is provided, it loads it into PyMOL,
        applies the mutation script, and saves the rendered PNG inside PyMOL_rendering/.
        """
        script_content = self.generate_pymol_script(annotated_variants)
        script_path = os.path.join(self.output_dir, "render_mutations.pml")

        # Save the script to disk
        with open(script_path, "w") as f:
            f.write(script_content)

        print(f"PyMOL script successfully saved to: {script_path}")

        try:
            import pymol
            # Check if pymol arguments are initialized to prevent second-launch deadlocks
            if not pymol:
                pymol.finish_launching(['pymol', '-qc'])

            pymol.cmd.delete("all")
            if PDB_file_path and os.path.exists(PDB_file_path):
                pymol.cmd.load(PDB_file_path, self.protein_object_name)
            else:
                pymol.cmd.fetch("1crn", self.protein_object_name, async_=0)

            pymol.cmd.do(script_content)

            image_output = os.path.abspath(os.path.join(self.output_dir, "mutation_render.png"))
            return f"Render successfully generated at: {image_output}"
        except ImportError:
            return f"PyMOL Python library not detected in environment."

    def render_rotation_video(self, annotated_variants: list, PDB_file_path: str = None, total_frames: int = 90,
                              fps: int = 24) -> str:
        """
        Renders a 360-degree rotation frame sequence via PyMOL in 16:9 landscape and compiles them into an MP4 video.
        """
        frames_dir = os.path.abspath(os.path.join(self.output_dir, "rotation_frames"))
        os.makedirs(frames_dir, exist_ok=True)

        output_video_path = os.path.abspath(os.path.join(self.output_dir, "mutation_rotation_360.mp4"))

        affected_residues = set()
        for v in annotated_variants:
            pos = v["position"]
            aa_pos = (pos // 3) + 1
            affected_residues.add(str(aa_pos))
        resi_string = "+".join(sorted(affected_residues, key=int)) if affected_residues else "1"

        try:
            import pymol
            if not pymol:
                pymol.finish_launching(['pymol', '-qc'])

            pymol.cmd.delete("all")
            if PDB_file_path and os.path.exists(PDB_file_path):
                pymol.cmd.load(PDB_file_path, self.protein_object_name)
            else:
                pymol.cmd.fetch("1crn", self.protein_object_name, async_=0)

            # Setup styling and focus
            pymol.cmd.hide('everything', self.protein_object_name)
            pymol.cmd.show('cartoon', self.protein_object_name)
            pymol.cmd.color('grey80', self.protein_object_name)
            pymol.cmd.select('mutated_sites', f'resi {resi_string} and {self.protein_object_name}')
            pymol.cmd.show('sticks', 'mutated_sites')
            pymol.cmd.color('firebrick', 'mutated_sites')

            # Set landscape 16:9 viewport resolution (1280x720)
            pymol.cmd.viewport(1280, 720)
            pymol.cmd.zoom('mutated_sites', buffer=5.0)

            angle_step = 360.0 / total_frames
            print(f"Rendering {total_frames} landscape rotation frames at {angle_step}° increments...")

            for i in range(total_frames):
                frame_path = os.path.join(frames_dir, f"frame_{i:03d}.png").replace("\\", "/")
                # Dump frames at 1280x720 landscape resolution
                pymol.cmd.png(frame_path, width=1280, height=720, ray=0)
                pymol.cmd.rotate('y', angle_step)

            # Safely quit PyMOL session at the very end of all rendering tasks
            pymol.cmd.quit()

            print("Compiling landscape frames into 360-degree MP4 video...")
            images = sorted([img for img in os.listdir(frames_dir) if img.endswith(".png")])
            if not images:
                return "Error: No frames found to compile into video."

            first_frame_path = os.path.join(frames_dir, images[0])
            sample_img = cv2.imread(first_frame_path)
            if sample_img is None:
                return f"Error: Could not read sample frame at {first_frame_path}"

            height, width, _ = sample_img.shape

            fourcc = cv2.VideoWriter_fourcc(*'mp4v')
            video_writer = cv2.VideoWriter(output_video_path, fourcc, fps, (width, height))

            for image in images:
                img_path = os.path.join(frames_dir, image)
                frame = cv2.imread(img_path)
                if frame is not None:
                    video_writer.write(frame)

            video_writer.release()
            return f"360-degree rotation video successfully compiled at: {output_video_path}"

        except ImportError as e:
            return f"Missing required library (PyMOL or OpenCV): {e}"