#!/pyvenv/bin/python

"""
Oxford Nanopore Technology (ONT) variant calling (vcall) utility script.

Ensure the directory 'raw_reads' containing fastq files present, otherwise this script will not run.

To setup the environment (directory creations and preparing reference) automatically
    ont-vcall setup

To run, per sample; use <sample_name> without the .fastq file extension
    ont-vcall run --sample <sample_name>

Please run on per gene segment basis!!
"""

import argparse
from subprocess import Popen, PIPE
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd


# ------------------------------------------------------------------------------
def prepare_environment() -> None:
    if Path("raw_reads").exists():
        print("  The 'raw_reads' directory is present, proceed...")
    else:
        print("  Could not find 'raw_reads' directory, aborting...")
        exit(1)

    for directory in ["output_results", "output_plots"]:
        if Path(directory).exists():
            print(f"  The folder '{directory}' already exists & I can't overwrite it, exiting...")
            exit(1)
        else:
            print(f"  Creating the folder '{directory}' for you...")
            Path.mkdir(directory)

    ACCESSION = "AF389118"

    print("")
    print("  Defaulting to 'PR8 HA (segment 4)'; GenBank AF389118")
    user_input_accession = input("  Enter alternative accession ID for reference (or leave blank): ")
    if len(user_input_accession) > 0:
        ACCESSION = user_input_accession

    for file in ["reference.fasta", "reference.gff"]:
        if Path(file).exists():
            print(f"  {file} exists, proceed...")
        else:
            if file == "reference.fasta":
                print(f"  '{file}' reference not found, pulling from GenBank...")
                pull_genbank(record="fasta", accession=ACCESSION)
            elif file == "reference.gff":
                print(f"  '{file}' feature file not found, pulling from GenBank...")
                pull_genbank(record="gff", accession=ACCESSION)


def pull_genbank(record: str, accession: str):
    if record == "fasta":
        command = f"bio fetch {accession} | bio fasta > reference.fasta"
        process = Popen(command, shell=True, stdout=PIPE)
        _, _ = process.communicate()
    elif record == "gff":
        command = f"bio fetch {accession} | bio gff > reference.gff"
        process = Popen(command, shell=True, stdout=PIPE)
        _, _ = process.communicate()
    else:
        print("  What??")
        exit(1)


def call_variants(sample: str) -> None:
    # This is a dirty way of doing it
    shell_script = Path("vcall_ont.sh")

    command_awk = "awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4}'"

    script_content = f"""
#!/usr/bin/env bash
set -xeu

amplicon-dist --fq raw_reads/{sample}.fastq
mv raw_reads/{sample}.fastq.jpg output_plots/{sample}_length_histogram.jpg

minimap2 -ax map-ont reference.fasta "raw_reads/{sample}.fastq" | samtools view --bam | samtools sort -o "{sample}.bam"

samtools index "{sample}.bam"

samtools mpileup --count-orphans --max-depth 600000 --no-BAQ --min-MQ 30 --fasta-ref reference.fasta "{sample}.bam" | ivar variants -p "{sample}" -q 20 -t 0.03 -r "reference.fasta" -g "reference.gff"

samtools mpileup {sample}.bam --fasta-ref reference.fasta | {command_awk} > output_plots/{sample}_coverage.tsv

mv {sample}.bam {sample}.bam.bai {sample}.tsv ./output_results
"""

    shell_script.write_text(script_content)
    process = Popen("bash vcall_ont.sh", shell=True, stdout=PIPE)
    _, _ = process.communicate()
    
    shell_script.unlink()


def plot_coverage(sample: str) -> None:
    _column_name = ["Chr", "Position", "Base", "Coverage"]
    df = pd.read_csv(f"output_plots/{sample}_coverage.tsv", sep="\t", header=None, names=_column_name)

    fig, ax = plt.subplots(figsize=(10, 3.5))

    max_y = df["Coverage"].max() * 1.25 + 25
    ax.plot(df["Position"], df["Coverage"], color="darkslategray")
    ax.fill_between(df["Position"], y1=df["Coverage"], y2=0, color="darkslategray", alpha=0.1)

    ax.set_ylim(0, max_y)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Coverage")
    ax.set_xlabel("Nucleotide position")
    ax.set_title(sample, loc="left", fontsize=10)

    for y in range(0, int(max_y), 50):
        ax.axhline(y=y, color="silver", linewidth=0.5, alpha=0.25, zorder=-100)

    fig.savefig(f"output_plots/{sample}_coverage.png", dpi=150, format="png", bbox_inches="tight")


# ------------------------------------------------------------------------------
def cli():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="subcommand")

    subparsers.add_parser("setup", help="Check if the environment is ready.")

    subparsers_run = subparsers.add_parser("run", help="Align and call variant.")
    subparsers_run.add_argument("--sample", "-s", type=str, required=True, help="Sample name without file extension (must exist in 'raw_reads' folder)")

    args = parser.parse_args()

    if args.subcommand == "setup":
        prepare_environment()
    elif args.subcommand == "run":
        call_variants(sample=args.sample)
        plot_coverage(sample=args.sample)
    else:
        parser.print_help()


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    cli()
