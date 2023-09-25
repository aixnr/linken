#!/bin/env bash

set -eu

bwa mem -M {{ .BwaIndex }} {{ .R1 }} {{ .R2 }} -t {{ .Threads }} | \
  samtools view --bam | \
  samtools sort -o {{ .Bwa }}

samtools index {{ .Bwa }}

