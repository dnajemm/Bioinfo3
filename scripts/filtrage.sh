#!/usr/bin/env bash
# Filters VCF files based on per-sample coverage and allele ratio
VCF_DIR="/data/projet3/Output/results_variant_calling"  # folder with input vcf files
OUT_DIR="/data/projet3/Output/results_filtering"    # folder for outputs



for vcf in "$VCF_DIR"/*.vcf.gz; do # Loop through each compressed VCF file
  base=$(basename "${vcf%.vcf.gz}")
  out="$OUT_DIR/${base}.filtered.vcf.gz"

  
  bcftools view -v snps -m2 -M2 "$vcf" \
  | bcftools filter -Oz -o "$out" -i '
      (FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 5 &&
      (FORMAT/AD[0:0] + FORMAT/AD[0:1]) > 0 &&
      (FORMAT/AD[0:1]) / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) >= 0.40 &&
      (FORMAT/AD[0:1]) / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) <= 0.60
    '

  bcftools index -f "$out"
done

