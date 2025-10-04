#!/usr/bin/env bash

REF="/data/projet3/reference/brbr.fasta"            # reference genome
Mapping_reads_DIR="/data/projet3/results_mapping"   # folder with bam files
OUT_DIR="/data/projet3/results_variant_calling"     # output folder
THREADS=8


for f in $Mapping_reads_DIR/*.bam; do
    SAMPLE=$(basename "$f" .bam)     # sample name
    echo "Processing $SAMPLE ..."
    
    bcftools mpileup -Ou -f $REF "$f" \
    | bcftools call -mv -Oz -o $OUT_DIR/${SAMPLE}.vcf.gz

    tabix -p vcf $OUT_DIR/${SAMPLE}.vcf.gz
done
