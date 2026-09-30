import os
import pymol
import numpy as np


def _kabsch_fit(mobile: np.ndarray, target: np.ndarray):
    """
    Compute the rotation R and translation t that best superpose `mobile`
    onto `target` (both Nx3 arrays, same length, index-paired) in a
    least-squares sense. Returns (R, t) such that `mobile @ R.T + t`
    approximates `target`.

    This is needed because two independently solved PDB structures (e.g.
    an apo and a holo form of the same protein) are essentially never in
    the same coordinate frame — each has its own arbitrary crystallographic
    origin/orientation. Interpolating between their raw coordinates without
    first superposing them produces a large rigid-body slide/rotation of
    the whole molecule instead of a clean local conformational change.
    """
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    mobile_c = mobile - mobile_center
    target_c = target - target_center

    H = mobile_c.T @ target_c
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = target_center - R @ mobile_center
    return R, t


def _apply_transform(coords: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    return coords @ R.T + t


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
                                       output_frames_dir: str, num_frames: int = 60,
                                       align_selection: str = "name CA") -> str:
        """
        Phases (matched to fraction of num_frames):
          0%   - 25%: ligand travels INTO the pocket, protein stays in the apo conformation
          25%  - 75%: ligand is bound/settled, protein conformation morphs apo -> holo -> apo
          75%  - 100%: ligand travels back OUT of the pocket, protein returns to apo
        """

        os.makedirs(output_frames_dir, exist_ok=True)
        pymol.cmd.delete("all")

        pymol.cmd.load(apo_protein, "apo_src")
        pymol.cmd.load(holo_protein, "holo_src")

        apo_records = []
        pymol.cmd.iterate_state(1, "apo_src", "apo_records.append((name, x, y, z))",
                                 space={'apo_records': apo_records})
        holo_records = []
        pymol.cmd.iterate_state(1, "holo_src", "holo_records.append((name, x, y, z))",
                                 space={'holo_records': holo_records})

        apo_names = [r[0] for r in apo_records]
        apo_coords = np.array([r[1:] for r in apo_records], dtype=float)
        holo_names = [r[0] for r in holo_records]
        holo_coords = np.array([r[1:] for r in holo_records], dtype=float)

        if len(apo_coords) != len(holo_coords):
            print(f"WARNING: apo ({len(apo_coords)} atoms) and holo ({len(holo_coords)} atoms) "
                  f"atom counts differ. Truncating to the shorter structure — if the two files "
                  f"aren't in identical atom order, the per-atom morph correspondence will be wrong. "
                  f"Consider renumbering/aligning the two PDBs' atom records before running this.")
            min_len = min(len(apo_coords), len(holo_coords))
            apo_names = apo_names[:min_len]
            holo_names = holo_names[:min_len]
            apo_coords = apo_coords[:min_len]
            holo_coords = holo_coords[:min_len]

        name_filter = align_selection.replace("name ", "").strip() if align_selection.startswith("name ") else None
        if name_filter:
            fit_mask = [i for i in range(len(apo_names))
                        if apo_names[i] == name_filter and holo_names[i] == name_filter]
        else:
            fit_mask = list(range(len(apo_names)))

        if len(fit_mask) < 3:
            print(f"WARNING: only {len(fit_mask)} atoms matched align_selection={align_selection!r}; "
                  f"falling back to fitting on all matched atoms instead.")
            fit_mask = list(range(len(apo_names)))

        R, t = _kabsch_fit(holo_coords[fit_mask], apo_coords[fit_mask])
        holo_coords = _apply_transform(holo_coords, R, t)

        pymol.cmd.delete("apo_src")

        pymol.cmd.load(apo_protein, "protein_morph", state=1)
        for f in range(2, num_frames + 1):
            pymol.cmd.create("protein_morph", "protein_morph", 1, f)

        has_ligand = False
        if ligand_path and os.path.exists(ligand_path):
            pymol.cmd.load(ligand_path, "ligand_template")
            has_ligand = True
        else:

            pymol.cmd.select("potential_ligand", "holo_src and hetatm and not solvent")
            if pymol.cmd.count_atoms("potential_ligand") > 0:
                pymol.cmd.create("ligand_template", "potential_ligand")
                has_ligand = True
            pymol.cmd.delete("potential_ligand")

        pymol.cmd.delete("holo_src")

        if has_ligand:

            raw_ligand_coords = []
            pymol.cmd.iterate_state(1, "ligand_template", "raw_ligand_coords.append((x, y, z))",
                                     space={'raw_ligand_coords': raw_ligand_coords})
            raw_ligand_coords = np.array(raw_ligand_coords, dtype=float)
            aligned_ligand_coords = _apply_transform(raw_ligand_coords, R, t)

            _lig_idx = [0]

            def _seed_ligand_coords(x, y, z):
                i = _lig_idx[0]
                if i < len(aligned_ligand_coords):
                    _lig_idx[0] += 1
                    return list(aligned_ligand_coords[i])
                return [x, y, z]

            pymol.cmd.alter_state(1, "ligand_template", "(x, y, z) = _seed_ligand_coords(x, y, z)",
                                  space={'_seed_ligand_coords': _seed_ligand_coords})

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

        pymol.cmd.show("cartoon", "protein_morph")
        pymol.cmd.color("slate", "protein_morph")
        pymol.cmd.bg_color("white")
        pymol.cmd.set("ray_shadows", 0)
        pymol.cmd.set("antialias", 2)

        print(f"Rendering phase-controlled binding simulation across {num_frames} frames with locked static camera...")

        entry_offset = np.array([0.0, 0.0, 20.0])

        pymol.cmd.frame(1)
        if has_ligand:
            zoom_sele = "(protein_morph within 5 of ligand)"
            if pymol.cmd.count_atoms(zoom_sele) == 0:
                print("WARNING: no protein_morph atoms found within 5 of ligand after alignment; "
                      "widening the zoom cutoff to 10 and re-checking.")
                zoom_sele = "(protein_morph within 10 of ligand)"
            pymol.cmd.zoom(zoom_sele, buffer=2.5)
        else:
            pymol.cmd.zoom("protein_morph", buffer=4.0)

        locked_view = pymol.cmd.get_view()

        for i in range(1, num_frames + 1):
            fraction = (i - 1) / max(1, num_frames - 1)
            pymol.cmd.frame(i)

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

                pymol.cmd.hide("sticks", "protein_morph")
                pymol.cmd.color("slate", "protein_morph")

                # Active site highlighting & hydrogen bonds specific to state i
                pymol.cmd.delete("h_bonds")
                pymol.cmd.select(
                    "active_site",
                    f"(protein_morph within 4.5 of (ligand and state {i})) and state {i}"
                )
                pymol.cmd.show("sticks", "active_site")
                pymol.cmd.color("yellow", "active_site and elem C")

                try:
                    pymol.cmd.distance("h_bonds", f"(ligand and state {i})", f"(active_site and state {i})", 3.8)
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