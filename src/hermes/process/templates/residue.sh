#!/bin/env bash

set -eu

# eevee's substitute supports identifying amino acid substitution from variant table
# currently, substitute only supports for single gene at once
# multi-gene support is planned
# defaulting gene to Segment04-HA_AF389118 (referring the reference fasta)
# with transcript region within position 33 - 1727
# required the folder 'output_variant' to be created

eevee substitute \
  --sample {{ .Sample }} \
  --gene Segment04-HA_AF389118:33-1727 \
  --ref {{ .BwaIndex }}.fasta \
  --vcf output_variant/{{ .Sample }}_variant_table.tsv \
  --csv_dir output_variant
