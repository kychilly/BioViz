from bioviz import (
    MutationAligner,
    VariantImpactAnalyzer,
    GenomeTrackRenderer,
    DomainLollipopAnalyzer,
    StructuralIntegrator
)


def test_pipeline():
    reference_genome = "ATGGCTACCGGATCGATCGATCGATCGATCGATCGATCG"
    target_genome = "ATGGCTATCGGATCGATCGATCGATCGTTCGATCGATCG"

    print("1. Running Alignment...")
    aligner = MutationAligner(reference_seq=reference_genome)
    variants = aligner.analyze_mutations(target_genome)

    print("2. Running Variant Impact & Codon Mapping...")
    impact_analyzer = VariantImpactAnalyzer(reference_seq=reference_genome, target_seq=target_genome)
    annotated_variants = impact_analyzer.analyze_impact(variants)

    print("3. Generating Genomic Track & Timeline...")
    track_renderer = GenomeTrackRenderer(sequence_length=len(reference_genome), window_size=10)
    track_data = track_renderer.compute_density_and_hotspots(annotated_variants)
    print(track_renderer.render_ascii_timeline(track_data))

    print("4. Mapping Mutations to Protein Domains & Lollipop Architecture...")
    protein_length = len(reference_genome) // 3

    domain_analyzer = DomainLollipopAnalyzer(protein_length=protein_length)
    domain_mapping = domain_analyzer.map_mutations_to_domains(annotated_variants)
    print(domain_analyzer.render_ascii_lollipop_map(domain_mapping))

    print("5. Executing Structural Integration & Local PyMOL Rendering...")
    structural_integrator = StructuralIntegrator(protein_object_name="target_protein")
    render_status = structural_integrator.render_structure_image(annotated_variants)

    print(render_status)

    # print("6. Protein 360 degree video rendering...")
    # video_status = structural_integrator.render_rotation_video(annotated_variants)
    # print(video_status)
    #
    # print("7. Compiling existing PyMOL frames into 360-degree video...")
    # # Initialize the video compiler targeting your PyMOL rendering directory
    # VideoCompiler(
    #     frames_dir="src/bioviz/mutations/PyMOL_rendering/rotation_frames",
    #     output_dir="src/bioviz/mutations/PyMOL_rendering"
    # )
    print("8. Calculating Effect & Severity Heatmap Matrix...")
    from bioviz import EffectSeverityHeatmap

    heatmap_analyzer = EffectSeverityHeatmap()
    severity_matrix = heatmap_analyzer.compute_severity_matrix(annotated_variants)
    print(heatmap_analyzer.render_ascii_heatmap(severity_matrix))


if __name__ == "__main__":
    test_pipeline()