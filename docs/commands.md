# Standalone Execution Commands

Under the hood, Linken runs commands listed in this document.


## Running on an HPC Cluster

A quick example here is loading `bwa` tool.

```bash
# To find
module avail bwa

# To load
module load bwa
```

Of note, the linux command `nproc` can return an integer of processor core counts.

```bash
nproc
>> 8
```


## How To Read This Guide

Anything within square bracket `[]` is provided by the user.
Anything within curly braces `{}` is system-related variable.


## Building Reference Genome

Required tools:

- `bwa` for generating an index of a reference sequence.
- `samtools` also for generating an index of a reference sequence.
- `picard` for generating *sequence dictionary* of a reference sequence.

Mapping (a.k.a. assembling) reads requires a reference genome, typically in the `.fasta` (or `.fa`) format.
Therefore, an *index* of the reference genome has to be generated first using `bwa`, `samtools`, and `picard`.
We will come back to these tools again when it comes to mapping them.
As a rule of thumb, the reference `.fasta` genome and the *index* files generated from this step should reside in its own separate folder.

```bash
bwa index [reference.fasta] -p [prefix]
```

This command outputs 5 files ending with the following file extensions: `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`.
The files `.amb` and `.ann` are plaintext, whereas `.bwt`, `.pac`, and `.sa` are binary.

- `[reference.fasta]` refers to the reference genome to align against.
- `-p [prefix]` sets the output name; e.g., `-p SARS2` outputs all the files as `SARS2.{amb,ann,bwt,pac,sa}` instead of using the input file name. This is useful especially if the rest of the pipeline is harcoded with a singular input reference.

```bash
samtools faidx [reference.fasta]
```

This commands outputs a single file, `[reference.fasta].fai`.

```bash
java -jar ${PICARD_HOME}/picard.jar CreateSequenceDictionary \
    R=[reference.fasta] \
    O=[reference.dict]
```

This command outputs `[reference.dict]`.
On the NIAID HPC, substitute `${PICARD_HOME}` with `${EBROOTPICARD}`.

Information regarding the specific files generated from this step:

| File ext.       | Type      | From             | Description
| ---------       | ----      | ----             | :----------
| `.fasta`, `.fa` | plaintext | --               | Reference sequence obtained through database, e.g. NCBI GenBank.
| `.amb`          | plaintext | `bwa index`      |
| `.ann`          | plaintext | `bwa index`      |
| `.bwt`          | binary    | `bwa index`      |
| `.pac`          | binary    | `bwa index`      |
| `.sa`           | binary    | `bwa index`      |
| `.fasta.fai`    | plaintext | `samtools faidx` |
| `.dict`         | plaintext | `picard CreateSequenceDictionary` |

As a rule of thumb, use the same file prefix when generating the gnome reference.
I recommend using `ref.*` or `reference.*`.
The prefix is critical when performing mapping with `bwa mem` command.


## Preparing Data Generated by Sequencing Instrument

Before analyses can be done, reads from sequencing instruments have to be mapped to the reference genome.
In this step, we will start working with raw data generated from sequencing instrument.
The raw data is usually in the form of `*.fastq.gz` data.
The `.fastq` suffix means it is in `.fasta` format with annotated quality score, and the `.gz` suffix means it is compressed to save disk space.
Most bioinformatic tools can read straight from `*.gz` file and perform decompression at runtime.

Required tools:

- `fastqc` for generating quality check report. This tool generates a `.html` report which shows the features of the run. It also provides information on the adapter that's present on the reads.
- `cutadapt` for performing adapter trimming. Also includes quality control where it can remove low-quality bases.
- `bwa` for mapping reads back to the reference genome.
- `samtools` for converting the plaintext `.sam`-formatted output from `bwa` into `.bam` files (binary format) sort the alignments afterwards.

```bash
fastqc --threads [N] --nogroup [read.fastq.gz]
```

This command outputs `*_fastqc.html` report file and `*_fastqc.zip` scratch report files from `[read.fastq.gz]`.
No modification is done on the `[read.fastq.gz]` input file.
The parameter `--threads N` tells `fastqc` how many *threads* (i.e., cores) to use.
The `--nogroup` parameter disables grouping of bases for reads >50 bp.
During testing, running with `--nogroup` in user node caused `fastqc` to run into Java out-of-memory error right after the *analysis complete* step.
For more information about `fastqc` tool, run `fastqc -h` to show help.

```bash
cutadapt \
    --cores [N] \
    --error-rate 0.1 \
    --quality-cutoff 30 \
    -a [ADAPTER_3PRIME_FIRST_READ] -A [ADAPTER_3PRIME_SECOND_READ] \
    --output [R1.cutadapt.fastq] \
    --paired-output [R2.cutadapt.fastq] \
    [R1.fastq.gz] [R2.fastq.gz]

pigz --processes [N] --blocksize 2048 [R1.cutadapt.fastq] [R1.cutadapt.fastq.gz]
pigz --processes [N] --blocksize 2048 [R2.cutadapt.fastq] [R2.cutadapt.fastq.gz]
```

This command outputs `[R1.cutadapt.fastq]` and `[R2.cutadapt.fastq]` from `[R1.fastq.gz]` and `[R2.fastq.gz]`.
Run `cutadapt -h` to show help.

- `--cores N`: Specify the number of cores to run `cutadapt`.
- `--error-rate 0.1`:
- `--quality-cutoff 30`:
- `-a [ADAPTER_3PRIME_FIRST_READ]` and `-A [ADAPTER_3PRIME_SECOND_READ]`: The adapter sequence, depends on the kit. For example, Nextera kit would use `CTGTCTCTTATA` for both R1 and R2 (paired).
- `--output [R1.fastq]` and `--paired-output [R2.fastq]`

`cutadapt` does not compress its own `[R{n}.cutadapt.fastq]` output, so it is a good idea to compress it into `.gz` with `pigz`.
By default, `pigz` deletes the original files after compressing.
This behavior can be turned now by supplying `--keep` (or `-k`) to `pigz`.
Alternatively, the `fastqc` and `cutadapt` step can be replaced with `trimgalore`, see this section `HERE_SECTION_LINK`.

Once adapter sequence is removed, mapping can now be performed.

```bash
# Step 1 to generate [sample.bwa.bam] BAM file.
bwa mem -M -t [N] [path/to/reference/genome/][prefix] [R1.fastq.gz] [R2.fastq.gz] | \
    samtools view --threads [N] --bam | \
    samtools sort --threads [N] -o [sample.bwa.bam]

# Step 2 to generate [sample.bwa.bam].bai BAM index file.
samtools index --threads [N] [sample.bwa.bam]
```

This command outputs `[sample.bwa.bam]` and `[sample.bwa.bam].bai` files from `[R1.cutadapt.fastq.gz]` and `[R2.cutadapt.fastq.gz]`.

- `bwa mem -M -t N` performs alignment of the reads to the pre-generated reference sequence; `-M` to mark shorter split hits as secondary (for Picard compatibility) and `-t N` is the number of threads. Of note, `bwa mem` outputs plaintext `.sam`-formatted file to `stdout`. Here, we pipe `bwa mem` output to `samtools`.
- `samtools view --threads N --bam` converts from plaintext `.sam` format to binary `.bam` format as specified by the `--bam` parameter. Of note, `samtools view` also defaults to `stdout` output, which can be piped into `samtools sort`.
- `samtools sort --threads N -o [sample.bwa.bam]` sorts alignments by leftmost coordinates.
- `samtools index [sample.bwa.bam]` index the coordinate-sorted `[sample.bwa.bam]` file (for fast random access) and outputs `[sample.bwa.bam].bai`.
- All `samtools` subcommands used here can receive `--threads N` parameter to specify the number of threads to use.

At this point with both `sample.bwa.bam` and `[sample.bwa.bam].bai` present, the data can be analyzed with applications such as [Integrative Genome Viewer](https://software.broadinstitute.org/software/igv/) (maintained by the Broad Institute) or [Genome Workbench](https://www.ncbi.nlm.nih.gov/tools/gbench/) (maintained by the National Library of Medicine).
Note that both `sample.bwa.bam` and `[sample.bwa.bam].bai` files are needed, along with the reference sequence `[reference.fasta]` file that was used to generate the reference genome.


## Variant Calling

Required tools:

- `lofreq` for fast and sensitive variant-calling. The input is `[sample.bwa.bam]` and `[sample.bwa.bam].bai` along with `[reference.fasta]`, outputting `.vcf` variant call format file. The sub-command `call-parallel` has an implicit dependency, which is the output of `samtools index`.
- `bcftools` for filtering the called variants based on fixed threshold, e.g., `AF<0.005`, where `AF` stands for allelic frequency. By default, it outputs to `stdout`.
- `picard.jar` for recalibrating the base quality. The flag `-Xmx2g` tells `java` to allocate 2G as the maximum heap size.
- `GATK` for `HERE_WHY`.

```bash
# Call variants using lofreq
lofreq call-parallel --pp-threads [N] \
    -f [reference.fasta] \
    -o [variants.vcf] [sample.bwa.bam]

# Filter variants with allelic frequency above 0.005 with bcftools
# 2> to redirect stderr to /dev/null, 1> to redirect stdout to file
bcftools filter --exclude "AF<0.005" [variants.vcf] 2> /dev/null 1> [variants_filtered.vcf] 

# Recalibrate base qualities with picard
java -Xmx2g -jar $(PICARD) AddOrReplaceReadGroups \
    -I [sample.bwa.bam] \
    -O [sample_rg.bwa.bam] \
    -RGLB A -RGPL illumina -RGPU run -RGSM [sample]

# Recalibrate base with GATK
# First, generate .idx file with IndexFeatureFile
gatk IndexFeatureFile -I variants_filtered.vcf

gatk BaseRecalibrator \
    -R [reference.fasta] \
    -I [sample_rg.bwa.bam] \
    --known-sites [variants_filtered.vcf] \
    -O [recal_data.table]

# Generate new .bam
gatk ApplyBQSR \
  -R [reference.fasta] \
  -I [sample_rg.bwa.bam] \
  --bqsr-recal-file [recal_data.table] \
  -O [sample_rg_bqsr.bam]

# Final variant calling
lofreq call-parallel --pp-threads [N] \
    -f [reference.fasta] \
    -o [variants_recal.vcf] \
    [sample_rg_bqsr.bam]
```

Open the `[variants_recal.vcf]` file and pay a close attention to all of the variants that were called.


## Trim-Galore to Combine fastqc and cutadapt Calls

[Trim-Galore](https://github.com/FelixKrueger/TrimGalore) (GitHub: `FelixKrueger/TrimGalore`) is a wrapper for both `fastqc` and `cutadapt`, thus both tools must be available in `$PATH`.
Check whether the executable is available either as `trim_galore` or `trimgalore`.

```bash
trimgalore --cores [N] --paired R1.fastq.gz R2.fastq.gz
```

This command outputs `_val_1.fq.gz`, `.fastq.gz_trimming_report.txt`, `_val_2.fq.gz`, and `.gz_trimming_report.txt` from `R1.fastq.gz` and `R2.fastq.gz`.
If desired, the removal of the adapter can be confirmed by running `fastqc` tool on the `fq.gz` file generated by `trimgalore`.
In our case during testing, `trimgalore` detected the right adapter, which was `CTGTCTCTTATA` Nextera Transposase sequence.

- `--cores N`:
- `--paired`:

The nice thing about `trimgalore` is that it combines 3 calls into 1: `fastqc`, `cutadapt`, and `pigz` for quality check, trimming, and compressing.
Consider using `--output_dir` flag to set output directory.


## Checking for Read Coverage

IGV and Genome Workbench can show read coverage for your target genes.
Here is another approach using scripts.

```bash
samtools mpileup [sample.bwa.bam] --fasta-ref [reference.fasta] | \
    awk '{print $1"\t"$2"\t"$3"\t"$4}' > [sample]_coverage.tsv
```

The `mpileup` subcommand produces text pileup from an alignment.
The output is then piped into `awk` to filter out the relevant fields, where fields `$1`, `$2`, `$3`, and `$4` correspond to:

- `$1`: chromosome name
- `$2`: 1-based position on the chromosome
- `$3`: reference base from `reference.fasta`
- `$4`: number of reads at this position

