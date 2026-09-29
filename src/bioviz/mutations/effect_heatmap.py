import os
import csv
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


class EffectSeverityHeatmap:
    """Computes, visualizes, and exports a graphical image heatmap of mutation severity metrics."""

    def __init__(self, output_dir: str = "src/bioviz/mutations/heatmaps"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def compute_severity_matrix(self, annotated_variants: list) -> list:
        """
        Calculates or maps downstream impact metrics
        (Conservation score, Pathogenicity index, and Delta G stability change) for each variant.
        """
        matrix_rows = []

        for v in annotated_variants:
            pos = v.get("position", 0)
            ref = v.get("ref", "-")
            alt = v.get("alt", "-")
            aa_pos = (pos // 3) + 1

            # Simulated biological scores (hooks into SIFT, PolyPhen-2, FoldX, or ConSurf in production)
            conservation_score = round(0.75 + (pos % 23) * 0.01, 2)
            pathogenicity_index = round(0.45 + (pos % 31) * 0.015, 2)
            delta_g = round(-1.2 + (pos % 17) * 0.25, 2)

            # Categorize severity risk level
            if pathogenicity_index > 0.8 or abs(delta_g) > 2.0:
                severity = "HIGH"
            elif pathogenicity_index > 0.6 or abs(delta_g) > 1.0:
                severity = "MODERATE"
            else:
                severity = "BENIGN"

            matrix_rows.append({
                "residue_position": aa_pos,
                "mutation": f"{ref}{aa_pos}{alt}",
                "conservation": conservation_score,
                "pathogenicity": pathogenicity_index,
                "delta_g": delta_g,
                "severity": severity
            })

        # Automatically persist the data tables and generate the visual image heatmap
        self._export_matrix_to_files(matrix_rows)
        return matrix_rows

    def _export_matrix_to_files(self, matrix_rows: list):
        """Saves the calculated heatmap matrix as a CSV, text report, and visual PNG image."""
        if not matrix_rows:
            return

        # 1. Export as a structured CSV file
        csv_path = os.path.join(self.output_dir, "severity_matrix.csv")
        with open(csv_path, mode="w", newline="", encoding="utf-8") as csv_file:
            fieldnames = ["residue_position", "mutation", "conservation", "pathogenicity", "delta_g", "severity"]
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()
            for row in matrix_rows:
                writer.writerow(row)

        # 2. Export as a plain text report document
        txt_path = os.path.join(self.output_dir, "severity_matrix_report.txt")
        with open(txt_path, mode="w", encoding="utf-8") as txt_file:
            txt_file.write(self.render_ascii_heatmap(matrix_rows, strip_colors=True))

        # 3. Generate and export a true graphical image heatmap (.png)
        img_path = os.path.join(self.output_dir, "severity_heatmap.png")
        self._generate_image_heatmap(matrix_rows, img_path)

        print(f"Heatmap outputs generated successfully in '{self.output_dir}/':\n"
              f"  - Graphical Image: {os.path.abspath(img_path)}\n"
              f"  - CSV Data: {os.path.abspath(csv_path)}\n"
              f"  - Text Report: {os.path.abspath(txt_path)}")

    def _generate_image_heatmap(self, matrix_rows: list, output_path: str):
        """Uses Seaborn and Matplotlib to plot a clean, color-mapped grid image of the mutations."""
        plt.figure(figsize=(10, max(4, len(matrix_rows) * 0.6 + 1.5)))

        # Extract features for matrix plotting
        labels = [row["mutation"] for row in matrix_rows]
        data_matrix = np.array([
            [row["conservation"], row["pathogenicity"], row["delta_g"]]
            for row in matrix_rows
        ])

        # Plot heatmap using seaborn
        ax = sns.heatmap(
            data_matrix,
            annot=True,
            fmt=".2f",
            cmap="vlag",  # Diverging colormap (great for stability / changes)
            cbar_kws={'label': 'Metric Scale Value'},
            xticklabels=["Conservation", "Pathogenicity", "Delta G ($\Delta\Delta$G)"],
            yticklabels=labels,
            linewidths=0.5,
            linecolor="gray"
        )

        plt.title("Variant Effect & Severity Impact Heatmap", fontsize=14, fontweight="bold", pad=15)
        plt.xlabel("Scoring Metrics", fontsize=11, labelpad=10)
        plt.ylabel("Mutated Residues", fontsize=11, labelpad=10)
        plt.xticks(rotation=0)
        plt.yticks(rotation=0)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()

    def render_ascii_heatmap(self, matrix_rows: list, strip_colors: bool = False) -> str:
        """Generates a fallback text-based matrix report table."""
        if not matrix_rows:
            return "No variant impact data available to generate heatmap."

        lines = [
            "=" * 85,
            f"{'EFFECT & SEVERITY IMPACT HEATMAP MATRIX':^85}",
            "=" * 85,
            f"{'Residue':<10} | {'Mutation':<10} | {'Conservation':<14} | {'Pathogenicity':<14} | {'Delta G (kcal/mol)':<18} | {'Severity':<10}",
            "-" * 85
        ]

        COLORS = {} if strip_colors else {
            "HIGH": "\033[91m",
            "MODERATE": "\033[93m",
            "BENIGN": "\033[92m",
            "RESET": "\033[0m"
        }

        for row in matrix_rows:
            sev = row["severity"]
            color = COLORS.get(sev, "")
            reset = COLORS.get("RESET", "") if color else ""

            line = (
                f"{row['residue_position']:<10} | "
                f"{row['mutation']:<10} | "
                f"{row['conservation']:<14.2f} | "
                f"{row['pathogenicity']:<14.2f} | "
                f"{row['delta_g']:<18.2f} | "
                f"{color}{sev:<10}{reset}"
            )
            lines.append(line)

        lines.append("=" * 85)
        return "\n".join(lines)