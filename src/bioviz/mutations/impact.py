from Bio.Seq import Seq


class VariantImpactAnalyzer:
    """Maps raw variants to codons, amino acid substitutions, and functional severity categories."""

    def __init__(self, reference_seq: str, target_seq: str):
        self.reference_seq = str(reference_seq)
        self.target_seq = str(target_seq)

    def analyze_impact(self, variants: list) -> list:
        """
        Takes raw variant dictionaries, identifies affected codons, translates them,
        and classifies the functional consequence.
        """
        annotated_variants = []

        for v in variants:
            pos = v["position"] - 1  # Convert to 0-indexed position
            v_type = v["type"]

            impact_data = {
                "position": v["position"],
                "type": v_type,
                "reference_base": v["reference_base"],
                "mutated_base": v["mutated_base"],
                "codon_change": "N/A",
                "amino_acid_change": "N/A",
                "effect": "Unknown",
                "severity": "Low"
            }

            if v_type == "SNP":
                # Find the start of the 3-base codon block
                codon_start = (pos // 3) * 3

                if codon_start + 3 <= len(self.reference_seq) and codon_start + 3 <= len(self.target_seq):
                    ref_codon = self.reference_seq[codon_start:codon_start + 3]
                    mut_codon = self.target_seq[codon_start:codon_start + 3]

                    # Translate codons to amino acids using Biopython
                    ref_aa = str(Seq(ref_codon).translate())
                    mut_aa = str(Seq(mut_codon).translate())

                    impact_data["codon_change"] = f"{ref_codon} -> {mut_codon}"
                    impact_data["amino_acid_change"] = f"{ref_aa} -> {mut_aa}"

                    # Classify mutation effect
                    if ref_aa == mut_aa:
                        impact_data["effect"] = "Synonymous (Silent)"
                        impact_data["severity"] = "Low"
                    elif mut_aa == "*":
                        impact_data["effect"] = "Nonsense (Stop Gain)"
                        impact_data["severity"] = "High"
                    else:
                        impact_data["effect"] = "Missense"
                        impact_data["severity"] = "Moderate"

            elif v_type in ["Insertion", "Deletion"]:
                impact_data["effect"] = "Indel / Frameshift Risk"
                impact_data["severity"] = "High"

            annotated_variants.append(impact_data)

        return annotated_variants