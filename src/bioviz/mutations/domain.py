import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches


class DomainLollipopAnalyzer:
    """Maps mutations to functional protein domains, builds structural lollipop maps, and exports graphical image plots."""

    def __init__(self, protein_length: int, domains: list = None, output_dir: str = "src/bioviz/mutations/graphs"):
        self.protein_length = protein_length
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

        # Default functional domains mapped by amino acid residue coordinates
        self.domains = domains or [
            {"name": "N-Terminal Domain", "start": 1, "end": 5, "type": "Structural"},
            {"name": "Catalytic Active Site", "start": 6, "end": 10, "type": "Active Site"},
            {"name": "Transmembrane Region", "start": 11, "end": 20, "type": "Membrane"}
        ]

    def map_mutations_to_domains(self, annotated_variants: list) -> dict:
        """Maps annotated variants to predefined functional protein domains and triggers graphical export."""
        domain_report = {
            d["name"]: {
                "type": d["type"],
                "range": (d["start"], d["end"]),
                "mutations": [],
                "count": 0
            } for d in self.domains
        }
        unassigned_mutations = []

        for v in annotated_variants:
            pos = v["position"]
            aa_pos = (pos // 3) + 1  # Convert nucleotide coordinate to approximate amino acid position

            assigned = False
            for d_name, d_info in domain_report.items():
                start, end = d_info["range"]
                if start <= aa_pos <= end:
                    d_info["mutations"].append(v)
                    d_info["count"] += 1
                    assigned = True
                    break

            if not assigned:
                unassigned_mutations.append(v)

        mapping_result = {
            "domain_mapping": domain_report,
            "unassigned": unassigned_mutations
        }

        # Automatically export the graphical image plot and text report
        self._export_lollipop_artifacts(mapping_result)

        return mapping_result

    def _export_lollipop_artifacts(self, mapping_result: dict):
        """Exports the graphical PNG image and text report to the graphs directory."""
        img_path = os.path.join(self.output_dir, "lollipop_map.png")
        self._generate_image_lollipop(mapping_result, img_path)

        txt_path = os.path.join(self.output_dir, "lollipop_report.txt")
        with open(txt_path, "w", encoding="utf-8") as f:
            f.write(self.render_ascii_lollipop_map(mapping_result))

        print(f"Lollipop artifacts successfully exported to '{self.output_dir}/':\n"
              f"  - Graphical Image: {os.path.abspath(img_path)}\n"
              f"  - Text Report: {os.path.abspath(txt_path)}")

    def _generate_image_lollipop(self, mapping_result: dict, output_path: str):
        """Uses Matplotlib to draw a structural protein backbone, domain blocks, and mutation lollipops."""
        fig, ax = plt.subplots(figsize=(12, 5))

        # Draw protein backbone bar
        ax.plot([1, self.protein_length], [0, 0], color="black", linewidth=4, zorder=1)

        domain_colors = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3"]

        # Draw domain blocks along the backbone
        for i, (d_name, info) in enumerate(mapping_result["domain_mapping"].items()):
            start, end = info["range"]
            color = domain_colors[i % len(domain_colors)]

            rect = patches.Rectangle(
                (start, -0.15),
                max(1, end - start),
                0.3,
                facecolor=color,
                edgecolor="black",
                alpha=0.8,
                zorder=2
            )
            ax.add_patch(rect)

            # Add domain label text above the block
            ax.text(
                (start + end) / 2,
                0.35,
                f"{d_name}\n({info['type']})",
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
                color=color
            )

        # Draw lollipops for each assigned mutation
        for d_name, info in mapping_result["domain_mapping"].items():
            for v in info["mutations"]:
                pos = v.get("position", 0)
                aa_pos = (pos // 3) + 1
                ref = v.get("ref", "-")
                alt = v.get("alt", "-")
                mut_label = f"{ref}{aa_pos}{alt}"

                # Stem line
                ax.plot([aa_pos, aa_pos], [0, 1.0], color="firebrick", linestyle="--", linewidth=1.5, zorder=3)
                # Lollipop head marker
                ax.plot(aa_pos, 1.0, marker="o", markersize=9, color="firebrick", markeredgecolor="black", zorder=4)
                # Mutation text label
                ax.text(aa_pos, 1.15, mut_label, rotation=45, ha="left", va="bottom", fontsize=8, fontweight="bold")

        # Styling and limits
        ax.set_xlim(0, self.protein_length + 2)
        ax.set_ylim(-0.8, 2.2)
        ax.axhline(0, color='grey', linewidth=0.5)

        plt.title("Protein Domain Architecture & Mutation Lollipop Plot", fontsize=14, fontweight="bold", pad=25)
        plt.xlabel("Amino Acid Residue Position", fontsize=11, labelpad=10)
        ax.set_yticks([])  # Hide y-axis ticks

        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

    @staticmethod
    def render_ascii_lollipop_map(mapping_result: dict) -> str:
        """Renders an ASCII text-based structural lollipop map for local terminal viewing."""
        output = "\n--- Protein Domain Architecture & Lollipop Mapping ---\n"

        for domain_name, info in mapping_result["domain_mapping"].items():
            count = info["count"]
            d_type = info["type"]
            start, end = info["range"]

            # Generate lollipop symbols for each mutation hitting this domain
            lollies = " 🍭" * count if count > 0 else " [No Mutations]"
            output += f"Domain: {domain_name} [{start}-{end}] ({d_type}) ->{lollies} (Count: {count})\n"

        return output