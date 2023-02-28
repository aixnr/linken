# Standalone Execution Commands

Under the hood, Linken runs commands listed in this document.

## Testing the Pipeline

During testing, the pipeline was tested extensively against SARS-CoV-2 genome (~30 kb).
SARS-CoV-2 was chosen because it is a non-segmented virus, as opposed to influenza where it is a segmented virus.
A smaller non-segmented virus like VSV (11 kb) or any flavivirus could have been chosen for the testing, but using a virus with larger genomic size like SARS-CoV-2 would allow for the opportunity for measuring the performance of the pipeline.

## Running on an HPC Cluster

A quick example here is loading `bwa` tool.

```bash
# To find
module avail bwa

# To load
module load bwa
```

## Building Reference Genome

Required tools:

- `bwa` for `HERE_DESCRIPTION`.
- `samtools` for `HERE_DESCRIPTION`.
- `picard` for `HERE_DESCRIPTION`.

At the end of this step, the `HERE_FOLLOWING_FILES_ARE_GENERATED`.

```bash
bwa index [reference.fasta] -p [prefix]
```

This command outputs 5 files ending with the following file extensions: `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`.
The files `.amb` and `.ann` are plaintext, whereas `.bwt`, `.pac`, and `.sa` are binary.

- `[reference.fasta]` refers to the reference genome to align against.
- `-p [prefix]` sets the output name; e.g., `-p SARS2` outputs all the files as `SARS2.{amb,ann,bwt,pac,sa}` instead of using the input file name. This is useful especially if the rest of the pipeline is hardcoded to only use 1 input.

```bash
samtools faidx [reference.fasta]
```

This commands outputs a single file, `[reference.fasta].fai`.

```bash
java -jar ${PICARD_HOME}/picard.jar CreateSequenceDictionary \
    R=[reference.fasta] \
    O=[reference.dict]
```

On NIAID HPC, substitute `${PICARD_HOME}` with `${EBROOTPICARD}`.

Information regarding the specific files generated from this step:

| File ext.       | Type      | From             | Description
| ---------       | ----      | ----             | :----------
| `.fasta`, `.fa` | plaintext | --               |
| `.amb`          | plaintext | `bwa index`      |
| `.ann`          | plaintext | `bwa index`      |
| `.bwt`          | binary    | `bwa index`      |
| `.pac`          | binary    | `bwa index`      |
| `.sa`           | binary    | `bwa index`      |
| `.fasta.fai`    | plaintext | `samtools faidx` |
| `.dict`         | plaintext | `picard CreateSequenceDictionary` |
