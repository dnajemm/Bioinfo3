#!/bin/bash

VCF_DIR="/data/projet3/Output/results_variant_calling" # Directory containing input VCF files
OUT_DIR="/data/projet3/Output/results_filtering" # Directory for output filtered VCF files

mkdir -p "$OUT_DIR" #create output directory if it doesn't exist

for vcf in "$VCF_DIR"/*.vcf.gz; do
  base=$(basename "${vcf%.vcf.gz}") # Extract base name without extension
  out="$OUT_DIR/${base}.filtered.vcf.gz" # Define output file path and name

  # Filter for heterozygous SNPs with depth >= 5 and allele ratio between 0.4 and 0.6
  bcftools view -v snps -m2 -M2 "$vcf" \ 
  | bcftools filter -i '
      GT="het" &&
      (FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 5 &&
      (FORMAT/AD[0:1]) / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 0.4 &&
      (FORMAT/AD[0:1]) / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.6
    ' -Oz -o "$out"

  bcftools index -f "$out" # Index the filtered VCF file
done