# Linken NGS Pipeline for Influenza Deep Sequencing Analysis

Linken is derived from a previous influenza deep sequencing pipeline written by Matthew G. Angel at the Laboratory of Viral Diseases (LVD), National Institute of Allergy and Infectious Diseases (NIAID).
The original pipeline utilized common bioinformatic tools for preparing and for analyzing sequencing data generated from Illumina sequencing instruments.

Detailed documentations are available in the `docs` folder.

## Architecture Overview

## Consideration for Running on an HPC Cluster

It is likely a high-performance computing (HPC) does not support Docker or Podman for container orchestration.
The HPC available for NIAID investigators support [Singularity](https://sylabs.io/docs/), an open source container tool maintained by Sylabs.

`HERE_CONVERTING_TO_SINGULARITY_IMAGE`.

## Building A Reference Genome
