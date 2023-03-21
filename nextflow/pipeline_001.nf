ch_raw_reads = channel.fromFilePairs("raw_reads/*_{R1,R2}_*.fastq.gz")
bwa_index = "$projectDir/index/ref"

workflow {
  QC_REPORT(ch_raw_reads)
  TRIM(ch_raw_reads)
  MAP( TRIM.out )
  CALIBRATE_AND_CALL( MAP.out )
}


process QC_REPORT {
  publishDir "report/fastqc", mode: "copy"

  input:
    tuple val(sample_id), path(reads)
  
  output:
    path("*.html")

  script:
  """
  fastqc --threads \$(nproc) ${reads[0]} ${reads[1]}
  """
}


process TRIM {
  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("*.fq.gz")
  
  script:
  """
  trimgalore --paired ${reads[0]} ${reads[1]} --cores \$(nproc)
  """
}


process MAP {
  publishDir "data/alignment", mode: "copy"

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("*.bwa.bam")
    path("*.bai")

  script:
  """
  bwa mem -M ${bwa_index} ${reads[0]} ${reads[1]} -t \$(nproc) | \
    samtools view --bam | \
    samtools sort -o ${sample_id}.bwa.bam

  samtools index ${sample_id}.bwa.bam
  """
}


process CALIBRATE_AND_CALL {
  publishDir "data/variant", pattern: "*.vcf", mode: "copy"
  publishDir "data/alignment", pattern: "*.bam", mode: "copy"

  input:
    tuple val(sample_id), path(bwa)
    path(bai)

  output:
    path("*.vcf")
    path("*.bam")

  script:
  """
  lofreq call-parallel --pp-threads \$(nproc) \
    -f ${bwa_index}.fasta \
    -o ${sample_id}.vcf \
    ${bwa}

  bcftools filter \
    --exclude "AF<0.005" ${sample_id}.vcf \
    2> /dev/null 1> ${sample_id}_filtered.vcf

  java -Xmx2g -jar /usr/local/lib/picard.jar \
    AddOrReplaceReadGroups \
    -I ${bwa} -O ${sample_id}_rg.bwa.bam \
    -RGLB A -RGPL illumina -RGPU run -RGSM ${sample_id}

  gatk IndexFeatureFile -I ${sample_id}_filtered.vcf

  gatk BaseRecalibrator \
    -R ${bwa_index}.fasta \
    -I ${sample_id}_rg.bwa.bam \
    --known-sites ${sample_id}_filtered.vcf \
    -O ${sample_id}_recal-data.table

  gatk ApplyBQSR \
    -R ${bwa_index}.fasta \
    -I ${sample_id}_rg.bwa.bam \
    --bqsr-recal-file ${sample_id}_recal-data.table \
    -O ${sample_id}_rg_bqsr.bam

  lofreq call-parallel --pp-threads \$(nproc) \
    -f ${bwa_index}.fasta \
    -o ${sample_id}_filtered_recalibrated.vcf \
    ${sample_id}_rg_bqsr.bam
  """
}

