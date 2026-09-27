import os
import pymol
import numpy as np


class PyMoLEngine:
    def __init__(self, headless: bool = True):
        self.headless = headless
        # Safely launch PyMOL headless/GUI mode if not already running
        try:
            args = ['pymol', '-qc'] if self.headless else ['pymol']
            pymol.finish_launching(args)
        except Exception:
            pass  # Already launched or initialized in the current session

    def setup_scene(self, protein_path: str, ligand_path: str, output_frames_dir: str):
        os.makedirs(output_frames_dir, exist_ok=True)
        pymol.cmd.delete("all")

        pymol.cmd.load(protein_path, "protein")
        if ligand_path and os.path.exists(ligand_path):
            pymol.cmd.load(ligand_path, "ligand")
            pymol.cmd.show("sticks", "ligand")
            pymol.cmd.color("cyan", "ligand")

        pymol.cmd.show("cartoon", "protein")
        pymol.cmd.color("slate", "protein")
        pymol.cmd.bg_color("white")
        pymol.cmd.set("ray_shadows", 0)
        pymol.cmd.set("antialias", 2)
        pymol.cmd.zoom("protein")

    def render_full_binding_simulation(self, apo_protein: str, holo_protein: str, ligand_path: str,
                                       output_frames_dir: str, num_frames: int = 60) -> str:
        """
        Executes a phase-controlled binding trajectory across num_frames using robust
        frame-based state switching:
        - 0% to 25% (Phase 1): Ligand translates inward toward the active site.
        - 25% to 75% (Phase 2): Ligand resides in active site while protein morphs (Apo <-> Holo).
        - 75% to 100% (Phase 3): Ligand translates outward (unbinding).
        """
        os.makedirs(output_frames_dir, exist_ok=True)
        pymol.cmd.delete("all")

        # 1. Load aligned structures into multi-state object
        pymol.cmd.load(apo_protein, "protein_morph", state=1)
        pymol.cmd.load(holo_protein, "protein_morph", state=2)

        # Expand states to match num_frames if num_frames > 2
        if num_frames > 2:
            for i in range(3, num_frames + 1):
                pymol.cmd.create("protein_morph", "protein_morph", 2, i)

        # 2. Handle ligand loading
        has_ligand = False
        if ligand_path and os.path.exists(ligand_path):
            pymol.cmd.load(ligand_path, "ligand")
            has_ligand = True
        else:
            pymol.cmd.select("potential_ligand", "hetatm and not solvent")
            if pymol.cmd.count_selected("potential_ligand") > 0:
                pymol.cmd.create("ligand", "potential_ligand")
                has_ligand = True
            pymol.cmd.delete("potential_ligand")

        if has_ligand:
            pymol.cmd.show("sticks", "ligand")
            pymol.cmd.color("cyan", "ligand")

            base_coords = []
            pymol.cmd.iterate_state(1, "ligand", "base_coords.append((x, y, z))", space={'base_coords': base_coords})
            base_coords = np.array(base_coords)

        # General visual setup
        pymol.cmd.show("cartoon", "protein_morph")
        pymol.cmd.color("slate", "protein_morph")
        pymol.cmd.bg_color("white")
        pymol.cmd.set("ray_shadows", 0)
        pymol.cmd.set("antialias", 2)

        pymol.cmd.mset(f"1x{num_frames}")
        print(f"Rendering phase-controlled binding simulation across {num_frames} frames...")

        entry_offset = np.array([0.0, 0.0, 25.0])

        for i in range(1, num_frames + 1):
            fraction = (i - 1) / max(1, num_frames - 1)

            # Use PyMOL's frame command to drive both animation timeline and multi-state mapping
            pymol.cmd.frame(i)

            if has_ligand:
                if fraction <= 0.25:
                    t_factor = 1.0 - (fraction / 0.25)
                elif fraction <= 0.75:
                    t_factor = 0.0
                else:
                    t_factor = (fraction - 0.75) / 0.25

                current_shift = entry_offset * t_factor
                shifted_coords = base_coords + current_shift

                coord_idx = 0

                def update_coords(x, y, z):
                    nonlocal coord_idx
                    val = shifted_coords[coord_idx]
                    coord_idx += 1
                    return list(val)

                pymol.cmd.alter_state(
                    1,
                    "ligand",
                    "(x, y, z) = update_coords(x, y, z)",
                    space={'update_coords': update_coords, 'shifted_coords': shifted_coords}
                )

                pymol.cmd.delete("h_bonds")
                pymol.cmd.select("active_site", "protein_morph within 4.5 of ligand")
                pymol.cmd.show("sticks", "active_site")
                pymol.cmd.color("yellow", "active_site and elem C")

                try:
                    pymol.cmd.distance("h_bonds", "ligand", "active_site", 3.8)
                    pymol.cmd.color("hotpink", "h_bonds")
                except Exception:
                    pass

                zoom_target = "ligand or active_site"
            else:
                zoom_target = "protein_morph"

            pymol.cmd.zoom(zoom_target, buffer=4.0)
            rock_angle = 15.0 * np.sin(fraction * 2.0 * np.pi)
            if i == 1:
                pymol.cmd.orient(zoom_target)
            else:
                pymol.cmd.turn("y", rock_angle / num_frames)

            frame_path = os.path.join(output_frames_dir, f"frame_{i:04d}.png")
            pymol.cmd.png(frame_path, width=1200, height=900, ray=0)
            print(f"frame {i} completed")

        pymol.cmd.delete("all")
        print(f"All {num_frames} phase-controlled frames successfully exported.")
        return output_frames_dir