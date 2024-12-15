import sys
from pathlib import Path
from Bio import SeqIO
import tomllib
import csv
from substitution import Substitution


def config_read(tomlpath: str) -> dict:
    config_data: dict = {}

    if Path(tomlpath).exists():
        with open(tomlpath) as f:
            content = f.read()
            config_data = tomllib.loads(content)
    else:
        print(f"  The data file {tomlpath} not found...")
        print("  Exiting...")
        sys.exit()

    return config_data


def config_doctor(config_data: dict) -> None:
    csv_filename = config_data["sample_list"]
    reference_fasta = config_data["reference_fasta"]
    gene_name = config_data["gene_name"]
    tsv_vcf_directory = config_data["tsv_vcf_directory"]

    # Check if the 'paired.csv' (or user-specified sample file) exists
    if Path(csv_filename).exists():
        print(f"  The sample list '{csv_filename}' is present!")
    else:
        print(f"  The sample list '{csv_filename}' is absent!")
        print(f"  Aborting...")
        exit(1)

    # Check if the 'tsv_vcf_directory' is present or not
    if Path(tsv_vcf_directory).exists():
        print(f"  The '{tsv_vcf_directory}' directory is present!")
    else:
        print(f"  The '{tsv_vcf_directory}' directory is absent!")
        print(f"  Aborting...")
        exit(1)

    # Check if the reference fasta exist
    if Path(reference_fasta).exists():
        print(f"  The reference FASTA '{reference_fasta}' is present!")
    else:
        print(f"  The reference FASTA '{reference_fasta}' is absent!")
        print("  Aborting...")
        exit(1)

    # Check of the gene exists, requires biopython
    _record_id: list = []
    for _record in SeqIO.parse(reference_fasta, "fasta"):
        _record_id.append(_record.id)

    if gene_name in _record_id:
        print(f"  The gene name '{gene_name}' is identified!")
    else:
        print(f"  The gene name '{gene_name}' is not found!")
        print("  Aborting...")
        exit(1)

    # Print list of samples
    with open(csv_filename, "r") as f:
        print("  The following samples will be processed...")

        _csv_reader = csv.reader(f)
        for _row in _csv_reader:
            _sample_id = _row[0]
            _tsv_path = Path(f"{tsv_vcf_directory}/{_sample_id}_variant_table.tsv")

            if _tsv_path.exists():
                print(f"    {_sample_id}; {_tsv_path}")
            else:
                print(f"    {_sample_id}; could not locate its *_variant_table.tsv")


def substitution_bulk(config_data: dict) -> None:
    csv_filename = config_data["sample_list"]
    reference_fasta = config_data["reference_fasta"]
    gene_name = config_data["gene_name"]
    gene_info_str = str(config_data["cds"]).split("-")
    tsv_vcf_directory = config_data["tsv_vcf_directory"]

    kws_gene = {
        "gene_id": gene_name,
        "gene_info": [int(i) for i in gene_info_str]
    }

    delta = Substitution(reference_fasta_file=reference_fasta).tabularize(**kws_gene)

    with open(csv_filename, "r") as f:
        _csv_reader = csv.reader(f)
        for _row in _csv_reader:
            _sample_id = _row[0]
            _tsv_path = Path(f"{tsv_vcf_directory}/{_sample_id}_variant_table.tsv")
            print(f"  \nProcessing for {_tsv_path}")

            delta.map_variant(variant_table=_tsv_path).collect_residue_context()
            print(f"  Collected residue contexts for {_tsv_path}")

            delta.write_csv(f"{tsv_vcf_directory}/{_sample_id}_residue_table_{gene_name}.csv")
            print(f"  Wrote {_sample_id}_residue_table_{gene_name}.csv")

            delta.write_genbank(f"{tsv_vcf_directory}/{_sample_id}_variant_{gene_name}.gb")
            print(f"  Wrote {_sample_id}_variant_{gene_name}.gb")
