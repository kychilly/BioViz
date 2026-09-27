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
        Executes a phase-controlled binding trajectory with a rigidly locked static camera
        zoomed tightly onto the active site, clean single-conformation ligand rendering,
        and true protein structural morphing.
        """
        os.makedirs(output_frames_dir, exist_ok=True)
        pymol.cmd.delete("all")

        # 1. Load Apo and Holo structures for coordinate morphing
        pymol.cmd.load(apo_protein, "apo_src")
        pymol.cmd.load(holo_protein, "holo_src")

        apo_coords = []
        pymol.cmd.iterate_state(1, "apo_src", "apo_coords.append((x, y, z))", space={'apo_coords': apo_coords})
        apo_coords = np.array(apo_coords)

        holo_coords = []
        pymol.cmd.iterate_state(1, "holo_src", "holo_coords.append((x, y, z))", space={'holo_coords': holo_coords})
        holo_coords = np.array(holo_coords)

        pymol.cmd.delete("apo_src")
        pymol.cmd.delete("holo_src")

        if len(apo_coords) != len(holo_coords):
            min_len = min(len(apo_coords), len(holo_coords))
            apo_coords = apo_coords[:min_len]
            holo_coords = holo_coords[:min_len]

        # Create multi-state protein morph object
        pymol.cmd.load(apo_protein, "protein_morph", state=1)
        for f in range(2, num_frames + 1):
            pymol.cmd.create("protein_morph", "protein_morph", 1, f)

        # 2. Handle ligand loading cleanly
        has_ligand = False
        if ligand_path and os.path.exists(ligand_path):
            pymol.cmd.load(ligand_path, "ligand_template")
            has_ligand = True
        else:
            pymol.cmd.select("potential_ligand", "hetatm and not solvent")
            if pymol.cmd.count_selected("potential_ligand") > 0:
                pymol.cmd.create("ligand_template", "potential_ligand")
                has_ligand = True
            pymol.cmd.delete("potential_ligand")

        if has_ligand:
            # Create a multi-state ligand object matching num_frames to prevent artifacts
            pymol.cmd.create("ligand", "ligand_template", 1, 1)
            for f in range(2, num_frames + 1):
                pymol.cmd.create("ligand", "ligand_template", 1, f)
            pymol.cmd.delete("ligand_template")

            pymol.cmd.show("sticks", "ligand")
            pymol.cmd.show("spheres", "ligand")
            pymol.cmd.set("sphere_scale", 0.25, "ligand")
            pymol.cmd.color("cyan", "ligand")

            base_coords = []
            pymol.cmd.iterate_state(1, "ligand", "base_coords.append((x, y, z))", space={'base_coords': base_coords})
            base_coords = np.array(base_coords)

        # General visual styling
        pymol.cmd.show("cartoon", "protein_morph")
        pymol.cmd.color("slate", "protein_morph")
        pymol.cmd.bg_color("white")
        pymol.cmd.set("ray_shadows", 0)
        pymol.cmd.set("antialias", 2)

        print(f"Rendering phase-controlled binding simulation across {num_frames} frames with locked static camera...")

        entry_offset = np.array([0.0, 0.0, 20.0])

        # --- Establish 100% Static Camera Tightly Zoomed on Active Site ---
        pymol.cmd.frame(1)
        if has_ligand:
            # Center and zoom tightly on the active site pocket residues surrounding the docked ligand
            pymol.cmd.zoom("protein_morph within 5 of ligand", buffer=2.5)
        else:
            pymol.cmd.zoom("protein_morph", buffer=4.0)

        # Capture the immutable camera view matrix
        locked_view = pymol.cmd.get_view()

        for i in range(1, num_frames + 1):
            fraction = (i - 1) / max(1, num_frames - 1)
            pymol.cmd.frame(i)

            # --- Phase-Based Protein Structural Morphing (Apo -> Holo -> Apo) ---
            if fraction <= 0.25:
                morph_progress = 0.0
            elif fraction <= 0.75:
                morph_progress = (fraction - 0.25) / 0.50
            else:
                morph_progress = 1.0 - ((fraction - 0.75) / 0.25)

            interp_protein_coords = apo_coords * (1.0 - morph_progress) + holo_coords * morph_progress

            p_idx = 0

            def update_protein_coords(x, y, z):
                nonlocal p_idx
                if p_idx < len(interp_protein_coords):
                    val = interp_protein_coords[p_idx]
                    p_idx += 1
                    return list(val)
                return [x, y, z]

            pymol.cmd.alter_state(
                i,
                "protein_morph",
                "(x, y, z) = update_protein_coords(x, y, z)",
                space={'update_protein_coords': update_protein_coords}
            )

            if has_ligand:
                # --- Phase-Based Ligand Translation per state ---
                if fraction <= 0.25:
                    t_factor = 1.0 - (fraction / 0.25)
                elif fraction <= 0.75:
                    t_factor = 0.0
                else:
                    t_factor = (fraction - 0.75) / 0.25

                current_shift = entry_offset * t_factor
                shifted_coords = base_coords + current_shift

                l_idx = 0

                def update_ligand_coords(x, y, z):
                    nonlocal l_idx
                    if l_idx < len(shifted_coords):
                        val = shifted_coords[l_idx]
                        l_idx += 1
                        return list(val)
                    return [x, y, z]

                pymol.cmd.alter_state(
                    i,
                    "ligand",
                    "(x, y, z) = update_ligand_coords(x, y, z)",
                    space={'update_ligand_coords': update_ligand_coords}
                )

                # Active site highlighting & hydrogen bonds specific to state i
                pymol.cmd.delete("h_bonds")
                pymol.cmd.select("active_site", f"protein_morph and state {i} within 4.5 of ligand and state {i}")
                pymol.cmd.show("sticks", "active_site")
                pymol.cmd.color("yellow", "active_site and elem C")

                try:
                    pymol.cmd.distance("h_bonds", f"ligand and state {i}", f"active_site and state {i}", 3.8)
                    pymol.cmd.color("hotpink", "h_bonds")
                except Exception:
                    pass

            # STRICTLY ENFORCE LOCKED STATIC CAMERA ON EVERY FRAME
            pymol.cmd.set_view(locked_view)

            frame_path = os.path.join(output_frames_dir, f"frame_{i:04d}.png")
            pymol.cmd.png(frame_path, width=1200, height=900, ray=0)
            print(f"frame {i} completed")

        pymol.cmd.delete("all")
        print(f"All {num_frames} frames successfully exported with a locked active-site zoom.")
        return output_frames_dir