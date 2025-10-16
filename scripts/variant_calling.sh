#!/usr/bin/env bash
# Script for variant calling using bcftools: mpileup and bcftools: call on BAM files.
# takes as input the BAM files from the mapping step and a reference genome in fasta format 
# gives as output compressed VCF files containing the identified variants

REF="/data/projet3/Input/reference/brbr.fasta"            # reference genome
Mapping_reads_DIR="/data/projet3/Output/results_mapping"   # folder with bam files from mapping step
OUT_DIR="/data/projet3/Output/results_variant_calling"     # output folder
THREADS=8                                           # number of CPU cores used


for f in $Mapping_reads_DIR/*.bam; do  # loop through all .bam files
    SAMPLE=$(basename "$f" .bam)     # sample name
    
    
    bcftools mpileup -Ou -f $REF "$f" \   # generate genotype likelihoods -f reference genome 
    | bcftools call -mv -Oz -o $OUT_DIR/${SAMPLE}.vcf.gz 
      # call variants -mv : multiallelic and variant sites only -Oz : output in compressed VCF format

    tabix -p vcf $OUT_DIR/${SAMPLE}.vcf.gz
    # tabix index the VCF file for fast access
done
