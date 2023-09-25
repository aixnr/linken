#!/bin/env bash

set -eu

lofreq call-parallel --pp-threads {{ .Threads }} \
  -f {{ .BwaIndex }}.fasta \
  -o {{ .Sample }}.vcf \
  {{ .Bwa }}

bcftools filter \
  --exclude "AF<0.005" {{ .Sample }}.vcf \
  2> /dev/null 1> {{ .Sample }}_filtered.vcf

java -Xmx2g -jar /usr/local/lib/picard.jar \
  AddOrReplaceReadGroups \
  -I {{ .Bwa }} -O {{ .Sample }}_rg.bwa.bam \
  -RGLB A -RGPL illumina -RGPU run -RGSM {{ .Sample }}

gatk IndexFeatureFile -I {{ .Sample }}_filtered.vcf

gatk BaseRecalibrator \
  -R {{ .BwaIndex }}.fasta \
  -I {{ .Sample }}_rg.bwa.bam \
  --known-sites {{ .Sample }}_filtered.vcf \
  -O {{ .Sample }}_recal-data.table

gatk ApplyBQSR \
  -R {{ .BwaIndex }}.fasta \
  -I {{ .Sample }}_rg.bwa.bam \
  --bqsr-recal-file {{ .Sample }}_recal-data.table \
  -O {{ .Sample }}_rg_bqsr.bam

lofreq call-parallel --pp-threads {{ .Threads }} \
  -f {{ .BwaIndex }}.fasta \
  -o {{ .Sample }}_filtered_recalibrated.vcf \
  {{ .Sample }}_rg_bqsr.bam

# cleanup, move everything to output_map
mv {{ .Sample }}* output_map/

