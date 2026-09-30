import os
import pymol
import numpy as np
from bioviz.rendering.parser import MoleculeParser
from bioviz.rendering.compiler import VideoCompiler
from bioviz.rendering.engine import _kabsch_fit, _apply_transform


class ProteinMorphEngine:
    def __init__(self, headless: bool = True):
        self.headless = headless
        try:
            args = ['pymol', '-qc'] if self.headless else ['pymol']
            pymol.finish_launching(args)
        except Exception:
            pass

    def render_protein_morph(self, apo_protein: str, holo_protein: str,
                             output_frames_dir: str, num_frames: int = 60,
                             align_selection: str = "name CA"):

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
            fit_mask = list(range(len(apo_names)))

        R, t = _kabsch_fit(holo_coords[fit_mask], apo_coords[fit_mask])
        holo_coords = _apply_transform(holo_coords, R, t)

        pymol.cmd.delete("apo_src")
        pymol.cmd.delete("holo_src")

        pymol.cmd.load(apo_protein, "protein_morph", state=1)
        for f in range(2, num_frames + 1):
            pymol.cmd.create("protein_morph", "protein_morph", 1, f)

        pymol.cmd.show("cartoon", "protein_morph")
        pymol.cmd.color("slate", "protein_morph")
        pymol.cmd.bg_color("white")
        pymol.cmd.set("ray_shadows", 0)
        pymol.cmd.set("antialias", 2)

        pymol.cmd.frame(1)
        pymol.cmd.zoom("protein_morph", buffer=4.0)
        locked_view = pymol.cmd.get_view()

        print(f"Rendering protein conformational morph across {num_frames} frames...")

        for i in range(1, num_frames + 1):
            fraction = (i - 1) / max(1, num_frames - 1)
            pymol.cmd.frame(i)

            morph_progress = 0.5 * (1.0 - np.cos(fraction * 2 * np.pi))

            interp_protein_coords = apo_coords * (1.0 - morph_progress) + holo_coords * morph_progress

            p_idx = [0]

            def update_protein_coords(x, y, z):
                idx = p_idx[0]
                if idx < len(interp_protein_coords):
                    val = interp_protein_coords[idx]
                    p_idx[0] += 1
                    return list(val)
                return [x, y, z]

            pymol.cmd.alter_state(
                i,
                "protein_morph",
                "(x, y, z) = update_protein_coords(x, y, z)",
                space={'update_protein_coords': update_protein_coords}
            )

            pymol.cmd.set_view(locked_view)
            frame_path = os.path.join(output_frames_dir, f"frame_{i:04d}.png")
            pymol.cmd.png(frame_path, width=1200, height=900, ray=0)
            print(f"Rendered protein morph frame {i}/{num_frames}")

        pymol.cmd.delete("all")
        return output_frames_dir


def main():
    apo_id = "4ake"
    holo_id = "1ake"

    print(f"=== Step 1: Fetching Multi-State Structures ({apo_id.upper()} vs {holo_id.upper()}) ===")
    paths = MoleculeParser.fetch_multi_state_complex(apo_id=apo_id, holo_id=holo_id)

    aligned_apo = os.path.join("sample_data", "protein_apo_aligned.pdb")
    MoleculeParser.align_structures(
        mobile_pdb=paths["apo_protein"],
        target_pdb=paths["holo_protein"],
        output_aligned_path=aligned_apo
    )

    print("=== Step 2: Initializing Protein Morph Engine ===")
    output_frames_dir = "protein_morph_frames"
    engine = ProteinMorphEngine(headless=True)
    engine.render_protein_morph(
        apo_protein=aligned_apo,
        holo_protein=paths["holo_protein"],
        output_frames_dir=output_frames_dir,
        num_frames=20
    )

    print("=== Step 3: Compiling into rendering_videos/protein_morph.mp4 ===")
    output_dir = "rendering_videos"
    os.makedirs(output_dir, exist_ok=True)
    output_video = os.path.join(output_dir, "protein_morph.mp4")

    VideoCompiler.create_video_from_frames(output_frames_dir, output_video, fps=5)
    print(f"\n🎉 Protein morph video successfully generated at: {os.path.abspath(output_video)}")


if __name__ == "__main__":
    main()