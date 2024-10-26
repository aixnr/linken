#!/bin/env bash

set -eu

# samtools mpileup output
# $1 ==> chromosome name
# $2 ==> 1-based position on the chromosome
# $3 ==> reference sequence from {{ .BwaIndex }}.fasta (refseq)
# $4 ==> number of reads at this position

samtools mpileup output_map/{{ .Bwa }} \
  --fasta-ref {{ .BwaIndex }}.fasta | \
  awk '{print $1"\t"$2"\t"$3"\t"$4}' > output_coverage/{{ .Sample }}_coverage.tsv
