import os
from bioviz.mut_reading import (
    MutationParser, StructuralCoordinateMapper,
    StructureMutator, MutationImpactAnalyzer,
    MutationVisualizer, InteractionNetworkAnalyzer, MutationReporter
)


def test_mutation_pipeline():
    print("=== Testing Full Bioviz Mutation Pipeline (Steps 1-7) ===")

    # 1. Load mutation records
    records = MutationParser.load_mutations()
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Locate predetermined structure file (1crn.cif)
    structure_path = os.path.join(script_dir, "bioviz", "mut_reading", "predetermined_proteins", "1crn.cif")
    if not os.path.exists(structure_path):
        # Fallback relative check depending on execution directory
        structure_path = os.path.join(script_dir, "predetermined_proteins", "1crn.cif")

    if os.path.exists(structure_path):
        print(f"\n[INFO] Found structure file: {structure_path}")

        # 2. Map coordinates & residue offsets
        mapper = StructuralCoordinateMapper(structure_file=structure_path, chain_id="A", residue_offset=0)
        mapper.map_mutations(records)

        # Save baseline Wild-Type structure backup
        output_dir = os.path.join(script_dir, "output")
        wt_pdb_out = os.path.join(output_dir, "wt_1crn_clean.pdb")
        MutationReporter.save_structure(mapper.structure, wt_pdb_out)

        # 3. Perform in silico mutagenesis and rotamer optimization (Position 1: T -> A)
        target_pos = 1
        mutator = StructureMutator(structure=mapper.structure, chain_id="A")
        mut_res_result = mutator.mutate_residue(residue_position=target_pos, target_wt="T", mut_residue_1le="A")

        # Save Mutated Structure PDB
        mut_pdb_out = os.path.join(output_dir, f"mut_1crn_pos{target_pos}_A.pdb")
        MutationReporter.save_structure(mapper.structure, mut_pdb_out)

        # 4. Analyze structural and physicochemical impact
        analyzer = MutationImpactAnalyzer(structure=mapper.structure, chain_id="A")
        impact_result = analyzer.analyze_impact(residue_position=target_pos, wt_res="T", mut_res="A")

        # Merge results into a master record dictionary for reporting
        master_record = {**mut_res_result, **impact_result}

        # 5. Render PyMOL Comparative Image with Top-Left Color Key Overlay
        img_output = os.path.join(output_dir, f"mutation_{target_pos}_comparison.png")
        try:
            MutationVisualizer.render_mutation_comparison(
                wt_pdb_path=wt_pdb_out,
                mut_pdb_path=mut_pdb_out,
                residue_position=target_pos,
                chain_id="A",
                output_image_path=img_output
            )
        except Exception as e:
            print(f"[INFO] PyMOL rendering note or graphics unavailable: {e}")

        # 6. Compute Interaction Networks & Generate Fingerprint Difference Profile
        net_analyzer = InteractionNetworkAnalyzer(structure=mapper.structure, chain_id="A")
        wt_contacts = net_analyzer.compute_contacts(residue_position=target_pos)
        mut_contacts = net_analyzer.compute_contacts(residue_position=target_pos)

        plot_output = os.path.join(output_dir, f"interaction_fingerprint_{target_pos}.png")
        InteractionNetworkAnalyzer.generate_fingerprint_plot(
            wt_contacts=wt_contacts,
            mut_contacts=mut_contacts,
            residue_position=target_pos,
            output_path=plot_output
        )

        # 7. Compile and Export Comprehensive JSON & CSV Reports
        MutationReporter.compile_reports([master_record], output_dir=output_dir)

        print("\n[COMPLETE] Full pipeline execution completed successfully! All files generated in /output.")

    else:
        print(f"\n[ERROR] Structure file {structure_path} could not be found.")


if __name__ == "__main__":
    test_mutation_pipeline()