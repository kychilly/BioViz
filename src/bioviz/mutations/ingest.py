import os
from Bio import Entrez, SeqIO
import requests

# Set your contact email for NCBI API compliance
Entrez.email = "jyam478@gmail.com"


class GenomeIngestor:
    """Universal ingestion engine to pull or load genome sequences from any source."""

    @staticmethod
    def fetch_from_ncbi(accession_id: str, file_format: str = "fasta") -> str:
        """Programmatically fetches a sequence from NCBI Entrez using an accession ID."""
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype=file_format, retmode="text")
            seq_data = handle.read()
            handle.close()
            return seq_data
        except Exception as e:
            raise RuntimeError(f"Failed to fetch accession {accession_id} from NCBI: {e}")

    @staticmethod
    def fetch_from_ensembl(species: str, region: str) -> dict:
        """Programmatically fetches genomic region data from the Ensembl REST API."""
        url = f"https://rest.ensembl.org/sequence/region/{species}/{region}?content-type=application/json"
        response = requests.get(url, headers={"Content-Type": "application/json"})
        if response.status_code == 200:
            return response.json()
        else:
            raise RuntimeError(f"Ensembl API error: {response.status_code} - {response.text}")

    @staticmethod
    def load_local_file(file_path: str):
        """Parses local FASTA, FASTQ, or VCF files from disk."""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Local file not found: {file_path}")

        ext = file_path.split('.')[-1].lower()
        if ext in ["fasta", "fa"]:
            return list(SeqIO.parse(file_path, "fasta"))
        elif ext in ["fastq", "fq"]:
            return list(SeqIO.parse(file_path, "fastq"))
        elif ext == "vcf":
            with open(file_path, 'r') as f:
                return [line.strip() for line in f if not line.startswith('##')]
        else:
            raise ValueError(f"Unsupported file extension: .{ext}")