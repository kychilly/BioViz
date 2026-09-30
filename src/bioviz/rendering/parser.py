import os
import urllib.request
import pymol


class MoleculeParser:
    @staticmethod
    def validate_file_path(file_path: str) -> str:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Target PDB file not found at: {file_path}")
        return os.path.abspath(file_path)

    # Download PDB, clean polymer, extract organic ligands
    @staticmethod
    def process_arbitrary_pdb(pdb_id: str, output_dir: str = "sample_data") -> dict:
        os.makedirs(output_dir, exist_ok=True)
        raw_pdb = os.path.join(output_dir, f"{pdb_id.lower()}.pdb")
        protein_clean = os.path.join(output_dir, f"{pdb_id.lower()}_protein.pdb")
        ligand_path = os.path.join(output_dir, f"{pdb_id.lower()}_ligand.pdb")

        if not os.path.exists(raw_pdb):
            url = f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb"
            print(f"Fetching structure {pdb_id.upper()} from RCSB...")
            urllib.request.urlretrieve(url, raw_pdb)

        pymol.finish_launching(['pymol', '-qc'])
        pymol.cmd.delete("all")
        pymol.cmd.load(raw_pdb, "complex")

        pymol.cmd.remove("solvent")
        pymol.cmd.save(protein_clean, "complex and polymer")

        pymol.cmd.select("ligand_sel", "complex and not polymer and not solvent and not inorganic")
        count = pymol.cmd.count_atoms("ligand_sel")

        if count == 0:
            print(f"Warning: No organic ligand automatically detected in {pdb_id.upper()}.")
            ligand_path = None
        else:
            pymol.cmd.save(ligand_path, "ligand_sel")
            print(f"Successfully extracted ligand ({count} atoms) for {pdb_id.upper()}.")

        pymol.cmd.delete("all")

        return {
            "protein": protein_clean,
            "ligand": ligand_path,
            "raw_pdb": raw_pdb
        }

    @staticmethod
    def fetch_multi_state_complex(apo_id: str = "4ake", holo_id: str = "1ake", output_dir: str = "sample_data") -> dict:
        os.makedirs(output_dir, exist_ok=True)

        apo_raw = os.path.join(output_dir, f"{apo_id.lower()}_apo.pdb")
        holo_raw = os.path.join(output_dir, f"{holo_id.lower()}_holo.pdb")

        apo_clean = os.path.join(output_dir, "protein_apo_clean.pdb")
        holo_clean = os.path.join(output_dir, "protein_holo_clean.pdb")
        ligand_path = os.path.join(output_dir, "ligand_bound.sdf")

        # Download Apo State
        if not os.path.exists(apo_raw):
            url_apo = f"https://files.rcsb.org/download/{apo_id.lower()}.pdb"
            print(f"Fetching Apo state ({apo_id.upper()})...")
            urllib.request.urlretrieve(url_apo, apo_raw)

        # Download Holo State
        if not os.path.exists(holo_raw):
            url_holo = f"https://files.rcsb.org/download/{holo_id.lower()}.pdb"
            print(f"Fetching Holo state ({holo_id.upper()})...")
            urllib.request.urlretrieve(url_holo, holo_raw)

        # Process structures via PyMOL headless engine
        pymol.finish_launching(['pymol', '-qc'])
        pymol.cmd.delete("all")

        pymol.cmd.load(apo_raw, "apo")
        pymol.cmd.remove("solvent")
        pymol.cmd.save(apo_clean, "apo and polymer")

        pymol.cmd.load(holo_raw, "holo")
        pymol.cmd.remove("solvent")
        pymol.cmd.save(holo_clean, "holo and polymer")

        pymol.cmd.select("ligand_sel", "holo and not polymer and not solvent and not inorganic")
        pymol.cmd.save(ligand_path, "ligand_sel")

        pymol.cmd.delete("all")
        print("Successfully acquired and cleaned multi-state structural files.")

        return {
            "apo_protein": apo_clean,
            "holo_protein": holo_clean,
            "ligand": ligand_path
        }

    # Aligns protein backbone of PDB to target PDB
    @staticmethod
    def align_structures(mobile_pdb: str, target_pdb: str, output_aligned_path: str) -> str:

        pymol.finish_launching(['pymol', '-qc'])
        pymol.cmd.delete("all")

        pymol.cmd.load(target_pdb, "target_struct")
        pymol.cmd.load(mobile_pdb, "mobile_struct")

        alignment_result = pymol.cmd.super("mobile_struct and name CA", "target_struct and name CA")
        rmsd_score = alignment_result[0] if isinstance(alignment_result, tuple) else alignment_result

        print(f"Structural alignment complete. Backbone RMSD: {rmsd_score:.2f} Å")

        pymol.cmd.save(output_aligned_path, "mobile_struct")
        pymol.cmd.delete("all")

        return output_aligned_path