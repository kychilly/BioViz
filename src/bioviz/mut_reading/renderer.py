import os
import cv2
import numpy as np
from pymol import cmd


class MutationVisualizer:
    """Handles PyMOL-based visualization, superimposition, and ray-traced image rendering of mutations with post-processing legend overlay."""

    @staticmethod
    def render_mutation_comparison(
            wt_pdb_path: str,
            mut_pdb_path: str,
            residue_position: int,
            chain_id: str = "A",
            output_image_path: str = "output/mutation_site_comparison.png"
    ) -> str:
        """
        Loads wild-type and mutant structures, aligns them via superimposition,
        styles mutation sites with colored sticks, focuses the camera, exports a ray-traced image,
        and post-processes the output to overlay a clean color key in the top-left corner.
        """
        # Ensure output directory exists
        out_dir = os.path.dirname(output_image_path)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        # 1. Initialize PyMOL session safely
        cmd.reinitialize()

        # 2. Load Structures
        wt_obj = "wild_type"
        mut_obj = "mutant"
        cmd.load(wt_pdb_path, wt_obj)
        cmd.load(mut_pdb_path, mut_obj)

        # Retrieve residue selection reference
        mut_selection = f"resi {residue_position} and chain {chain_id}"

        # 3. Structural Superimposition (Align mutant backbone to wild-type backbone)
        cmd.super(f"{mut_obj} and chain {chain_id}", f"{wt_obj} and chain {chain_id}")

        # 4. Global Visual Styling
        cmd.hide("everything")
        cmd.show("cartoon", "all")

        # Color the rest of the protein neutral gray/slate tones
        cmd.color("gray80", wt_obj)
        cmd.color("slate", mut_obj)

        # 5. Style the Mutation Site as Distinct Colored Sticks
        cmd.show("sticks", f"{wt_obj} and {mut_selection}")
        cmd.show("sticks", f"{mut_obj} and {mut_selection}")

        # Color Wild-Type site cyan, Mutant site magenta
        cmd.color("cyan", f"{wt_obj} and {mut_selection}")
        cmd.color("magenta", f"{mut_obj} and {mut_selection}")

        # 6. Highlight interactions/distances at the mutation site
        cmd.distance(
            "mut_shift_vector",
            f"{wt_obj} and {mut_selection} and name CA",
            f"{mut_obj} and {mut_selection} and name CA"
        )
        cmd.color("yellow", "mut_shift_vector")

        # 7. Zoom tightly onto the mutated residue bounding box
        cmd.zoom(f"resi {residue_position} and chain {chain_id}", buffer=4.0)

        # 8. High-Resolution Ray Traced Image Export with Dark Publication Theme
        cmd.bg_color("black")
        cmd.set("ray_shadows", "on")
        cmd.set("antialias", 2)

        cmd.ray(2400, 1800)
        cmd.png(output_image_path)

        # 9. Post-Processing: Overlay Color Key Legend via OpenCV
        try:
            img = cv2.imread(output_image_path)
            if img is not None:
                h, w, _ = img.shape

                # Define legend box layout coordinates in the top-left corner
                box_x1, box_y1 = int(w * 0.02), int(h * 0.02)
                box_x2, box_y2 = int(w * 0.30), int(h * 0.24)

                # Draw solid dark container box with a crisp white border
                overlay = img.copy()
                cv2.rectangle(overlay, (box_x1, box_y1), (box_x2, box_y2), (15, 15, 15), -1)
                cv2.rectangle(overlay, (box_x1, box_y1), (box_x2, box_y2), (255, 255, 255), 3)

                # Blend for clean opacity
                alpha = 0.85
                cv2.addWeighted(overlay, alpha, img, 1 - alpha, 0, img)

                # Typography setup scaled to high resolution
                font = cv2.FONT_HERSHEY_DUPLEX
                font_scale = w / 1600.0
                thickness = max(2, int(w / 900))

                text_x = box_x1 + int(w * 0.018)
                curr_y = box_y1 + int(h * 0.04)
                line_spacing = int(h * 0.042)

                # Header: COLOR KEY
                cv2.putText(img, "COLOR KEY", (text_x, curr_y), font, font_scale * 1.1, (255, 255, 255), thickness + 1,
                            cv2.LINE_AA)
                curr_y += int(line_spacing * 1.1)

                # Divider rule
                cv2.line(img, (text_x, curr_y - int(line_spacing * 0.35)),
                         (box_x2 - int(w * 0.018), curr_y - int(line_spacing * 0.35)), (200, 200, 200), 2)

                # Key Items (BGR Color format)
                cv2.putText(img, "• Wild-Type (Cyan)", (text_x, curr_y), font, font_scale, (255, 255, 0), thickness,
                            cv2.LINE_AA)
                curr_y += line_spacing

                cv2.putText(img, "• Mutant (Magenta)", (text_x, curr_y), font, font_scale, (255, 0, 255), thickness,
                            cv2.LINE_AA)
                curr_y += line_spacing

                cv2.putText(img, "• Backbone (Cartoon)", (text_x, curr_y), font, font_scale, (200, 200, 200), thickness,
                            cv2.LINE_AA)
                curr_y += line_spacing

                cv2.putText(img, "• Clashes/Distance", (text_x, curr_y), font, font_scale, (0, 255, 255), thickness,
                            cv2.LINE_AA)

                # Save finalized image with embedded key
                cv2.imwrite(output_image_path, img)
                print(f"[SUCCESS] Color key successfully overlaid onto: {output_image_path}")
        except Exception as e:
            print(f"[WARNING] Post-processing legend overlay encountered an issue: {e}")

        print(f"[SUCCESS] High-resolution mutation comparison rendered and saved to: {output_image_path}")
        return output_image_path