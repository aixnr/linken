# Linken NGS Pipeline for Influenza Deep Sequencing Analysis

Linken is derived from a previous influenza deep sequencing pipeline written by Matthew G. Angel at the Laboratory of Viral Diseases (LVD), National Institute of Allergy and Infectious Diseases (NIAID).
The original pipeline utilized common bioinformatic tools for preparing and for analyzing sequencing data generated from Illumina sequencing instruments.

Detailed documentations are available in the `docs` folder.


## Architecture Overview

Designing a computational pipeline is inherently challenging.
A pipeline involves *moving* data from one bioinformatic program to another, and the sequence of action has to be right.
Making matter worse, sometimes different bioinformatic programs have different dependency requirements, resulting in conflict and program can crash.
To remedy the situation, the method that seems robust to build a stable bioinformatic computational pipeline is through containerization.
This often involves downloading the program's source codes, compiling in the target environment (e.g. `debian:11`), and convert the container into executable (i.e. Singular Image Format, SIF).
The resulting executable container is ready for usage.


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
The process is split into two parts:

1. Generating the index file for mapping with `bwa`, `samtools`, and `picard`.
2. Running `hermes` (defaults to Illumina with paired-end reads) for performing:
    - Quality control and base-trimming with `fastqc` and `trimgalore`.
    - Mapping and generating aligned reads with `bwa` and `samtools`.
    - Recalibrating reads with `picard` and `gatk`.
    - Variant calling with `lofreq`.

```bash
# Change directory (cd) into the analysis folder
cd analysis_folder

# Show linkenc lunar main interface
linkenc lunar

# Run lunar doctor to check for environment
linkenc lunar doctor

# List publicly available indexes (github.com/aixnr/linken-contrib; main branch)
linkenc lunar list

# Generate index files for a reference genome
linkenc lunar create --ref IAV-PR8-MtSinai
```

