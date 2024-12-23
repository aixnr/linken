# Linken NGS Pipeline for Viral Deep Sequencing Analysis

Linken is derived from a previous influenza deep sequencing pipeline written by Matthew G. Angel at the Laboratory of Viral Diseases (LVD), National Institute of Allergy and Infectious Diseases (NIAID).
The original pipeline utilized common bioinformatic tools for preparing and for analyzing sequencing data generated from Illumina sequencing instruments.

Detailed documentations are available in the `docs` folder.


## Architecture Overview

Designing a computational pipeline is inherently challenging.
A pipeline involves *moving* data from one bioinformatic program to another, and the sequence of action has to be right.
Making matter worse, sometimes different bioinformatic programs have different dependency requirements (*dependency hell*), resulting in conflict and program can crash.
*Containerization*, a strategy that is used here is one of several robust strategies that can be used to mitigate *dependency hell*.
However, Docker images are often unsupported on compute clusters.
Luckily, the Singularity/Apptainer project allows for *containers* to be safely executed on computer clusters.


## Consideration for Running on an HPC Cluster

It is likely a high-performance computing (HPC) does not support Docker or Podman for container orchestration.
The HPC available for NIAID investigators supports [Singularity](https://sylabs.io/docs/), an open source container tool by Sylabs.
Of note, Singular was renamed to [Apptainer](https://apptainer.org/) back in November 2021 (read: [community announcement](https://apptainer.org/news/community-announcement-20211130/)).
Apptainer is currently a project maintained by the [Linux Foundation](https://www.linuxfoundation.org/).

Singularity and Apptainer do support running containers generated by Docker, and basically containers generating from tools that adhere to the Open Containers Initiative (OCI) format.
Read and follow instructions on [Singularity](https://docs.sylabs.io/guides/2.6/user-guide/singularity_and_docker.html) and [Apptainer](https://apptainer.org/docs/user/main/docker_and_oci.html) documentation.

By default, Singularity/Apptainer mounts user's `$HOME` automatically inside the container, hence all files residing the user's `$HOME` is visible to the container.
This is not the case if the data is elsewhere, e.g., `/media/data` or `/hpcdata/`, i.e. anywhere outside of `$HOME`.
This will be a problem when running a Singularity/Apptainer container outside `$HOME` because the container won't be able to see files.
User can test this by issuing `echo $PWD`.

The fix is relatively simple.
Singularity honors environmental variable `SINGULARITY_BIND` (or `APPTAINER_BIND`).

```bash
# syntax for specific mount location
# specify source (src) and destination (dest)
export SINGULARITY_BIND=/src:/dest

# mount destination similar as source
export SINGULARITY_BIND=/hpcdata/
```

When `SINGULARITY_BIND` (e.g. to `/hpcdata/`), running `echo $PWD` from a container anywhere within `/hpcdata` will return the correct path.


## How To Use

Using Linken through Singularity/Apptainer container is the primary use case.
The process is split into three parts:

1. Generating the index file for mapping with `bwa`, `samtools`, and `picard`.
2. Running `hermes` (defaults to Illumina with paired-end reads) for performing:
    - Quality control and base-trimming with `fastqc` and `trimgalore`.
    - Mapping and generating aligned reads with `bwa` and `samtools`.
    - Recalibrating reads with `picard` and `gatk`.
    - Variant calling with `lofreq`.
3. Using `eevee` (which contains `coverage` and `substitution` analysis modules) to analyze the data.

```bash
# Change directory (cd) into the analysis folder
cd analysis_folder

# Show linken lunar main interface
linken lunar

# Run lunar doctor to check for environment
linken lunar doctor

# List publicly available indexes (github.com/aixnr/linken-contrib; main branch)
linken lunar list

# Generate index files for a reference genome
linken lunar create --ref IAV-PR8-MtSinai

# Process the raw_reads, assumes samples are defined in paired.csv by default
# The default behavior can be overridden
linken hermes
```

### Hermes

The default behavior for `hermes` when called without arguments is similar to this:

```bash
# hermes default behavior
hermes \
  -i paired.csv \
  -r 2 \
  -p illumina \
  -n 4
  -ref ./index/ref
```

Where

- `-i`: The input CSV specifying the sample names and their location. Defaults to `paired.csv`.
- `-r`: Legal input is either 1 or 2, corresponding to single or paired reads.
- `-p`: Sequencing platform. Only `illumina` is supported at the moment.
- `-n`: Number of processing threads to use, defaults to `4` threads.
- `-ref`: The reference files generated by `lunar create` command.


### Eevee

By default, only the `coverage` subcommand is called by `hermes` to generate plots showing the sequencing coverage and variant frequency, of which the outputs can be found at `output_plots` directory.
The `coverage` subcommand requires `sample_name` value, a `.tsv` generated by `samtools mpileup` which has the sequencing coverage information, and a `.tsv` of called variants (generated by `gatk VariantsToTable` from the `.vcf` file from `lofreq`):

```bash
linken eevee coverage \
  --sample {sample_name} \
  --cov_tsv output_coverage/{sample_name}_coverage.tsv \
  --vcf output_variant/{sample_name}_variant_table.tsv
```

`eevee` also supports *translating* the called nucleotide variants to amino acid through its `substitute` subcommand.
This portion of analysis needs to be specified by the user.

```bash
# signature
linken eevee substitute \
  --sample {sample_name} \
  --gene {gene_name:nn-mm} \
  --ref ./index/ref.fasta \
  --vcf output_variant/{sample_name}_variant_table.tsv \
  --csv_dir output_variant

# example
linken eevee substitute \
  --sample Sample01 \
  --gene Segment04-HA_AF389118:33-1727 \
  --ref ./index/ref.fasta \
  --vcf output_variant/Sample01_variant_table.tsv \
  --csv_dir output_variant
```

Crucially, the `--gene` parameter requires the gene name (e.g., `Segment04-HA_AF389118`) to precisely match with the gene present in the reference `.fasta`.
The number `:nn-mm` (e.g. `:33-1727`) refers to the open reading frame, of which the `substitute` subcommand will translate nucleotide codons to amino acid residues.

As of **2024/DEC/15**, `eevee` supports performing this analysis in bulk, provided that the `paired.csv` file (the same one used by `hermes`) is present.
This feature is delegated to the `substitute-bulk` subcommand, which takes either of these two parameters: `--check <configuration-file.toml>` or `--run <configuration-file.toml>`, where the user must supplied the `*.TOML` file containing the required parameters, shown below:

```toml
sample_list = "paired.csv"
reference_fasta = "index/ref.fasta"
tsv_vcf_directory = "output_variant"
gene_name = "Segment04-HA_AF389118"
cds = "33-1727"
```

To ensure the conditions are satisfied, issue the `check` mode, assuming the filename of the `*.TOML` file is `PR8.toml`

```bash
linken eevee substitute-bulk --check PR8.toml
```

If no error was returned, perform the analysis:

```bash
linken eevee substitute-bulk --run PR8.toml
```


### ont-vcall

As of `2024/DEC/22`, a utility script for performing variant calling on Oxford Nanopore Technology (ONT) sequencing reads was added; see `scripts/ont-vcall.py`.
This utility script relies on `minimap2` for read mapping and `iVar` for variant calling.

The hard requirement is the current directory must contain `raw_reads` sub-folder containing all the raw reads in `.fastq` format.

```bash
# setup the environment
# to pull reference sequence and feature file from GenBank
# and for folder creation ('output_results' and 'output_plots')
linken ont-vcall setup

# run, sample-by-sample, for 1 gene segment at a time
# the sample name must not contain the .fastq extension
# the sample must exist in the raw_reads sub-folder
linken ont-vcall run --sample <sample-name>
```
