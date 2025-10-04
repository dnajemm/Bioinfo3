#!/usr/bin/env bash
# Script for variant calling using bcftools

REF="/data/projet3/reference/brbr.fasta"            # reference genome
Mapping_reads_DIR="/data/projet3/results_mapping"   # folder with bam files
OUT_DIR="/data/projet3/results_variant_calling"     # output folder
THREADS=8


for f in $Mapping_reads_DIR/*.bam; do
    SAMPLE=$(basename "$f" .bam)     # sample name
    echo "Processing $SAMPLE "
    
    bcftools mpileup -Ou -f $REF "$f" \   # generate genotype likelihoods -f reference genome 
    | bcftools call -mv -Oz -o $OUT_DIR/${SAMPLE}.vcf.gz 
      # call variants -mv : multiallelic and variant sites only -Oz : output in compressed VCF format

    tabix -p vcf $OUT_DIR/${SAMPLE}.vcf.gz
    # tabix index the VCF file for fast access
done
