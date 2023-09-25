#!/bin/env bash

set -eu

# Using gatk VariantsToTable to turn a VCF into TSV containing called variants
# Only extracts for columns CHROM, POS, REF, and ALT
# And key AF in the INFO column into its own column

gatk VariantsToTable \
  --variant output_map/{{ .Sample }}_filtered_recalibrated.vcf \
  --fields CHROM \
  --fields POS \
  --fields REF \
  --fields ALT \
  --fields AF \
  --output output_variant/{{ .Sample }}_variant_table.tsv

# Generate plot using eevee's to show where the variants are in relation to the read coverage
# eevee coverage requires output_plots
# ./src/hermes/process/dir.go handles the generation of the output directory
eevee coverage \
  --sample {{ .Sample }} \
  --cov_tsv output_coverage/{{ .Sample }}_coverage.tsv \
  --vcf output_variant/{{ .Sample }}_variant_table.tsv

