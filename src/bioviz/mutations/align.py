from Bio import Align


class MutationAligner:
    """Performs sequence alignment and isolates variants against a reference sequence."""

    def __init__(self, reference_seq: str):
        self.reference_seq = str(reference_seq)
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -2
        self.aligner.extend_gap_score = -0.5

    def analyze_mutations(self, target_seq: str) -> list:
        """Aligns target sequence to reference and flags SNPs, insertions, and deletions."""
        target_seq = str(target_seq)
        alignments = self.aligner.align(self.reference_seq, target_seq)
        best_alignment = alignments[0]

        alignment_lines = format(best_alignment).split('\n')
        ref_aligned, target_aligned = alignment_lines[0], alignment_lines[2]

        mutations = []
        ref_pos = 0

        for idx, (r_char, t_char) in enumerate(zip(ref_aligned, target_aligned)):
            if r_char != '-':
                ref_pos += 1

            if r_char != t_char:
                mutation_type = "SNP" if r_char != '-' and t_char != '-' else (
                    "Insertion" if r_char == '-' else "Deletion")
                mutations.append({
                    "position": ref_pos,
                    "reference_base": r_char,
                    "mutated_base": t_char,
                    "type": mutation_type,
                    "alignment_index": idx
                })

        return mutations