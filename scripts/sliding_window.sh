#!/usr/bin/env bash
set -euo pipefail

IN_DIR="/data/projet3/Output/results_filtering"   # your *.filtered.vcf.gz files
OUT_DIR="/data/projet3/Output/results_1kb"        # output folder
mkdir -p "$OUT_DIR"

MERGED="$OUT_DIR/all_samples.het.1kb.snpden.tsv"   # one tidy file for R/py plots
echo -e "SAMPLE\tCHROM\tBIN_START\tN_HET\tSNP_DENSITY" > "$MERGED"

for vcf in "$IN_DIR"/*.filtered.vcf.gz; do
  s=$(basename "$vcf" .filtered.vcf.gz)           # e.g. YJS8112
  echo ">>> $s : 1 kb windows (het-only)"

  # het-only per sample (0/1 genotypes) → vcftools counts per-1kb bin
  # Use stdin pipe (vcftools --vcf -) to avoid the /dev/fd/63 error
  bcftools view -g het "$vcf" \
  | vcftools --vcf - \
             --SNPdensity 1000 \
             --out "$OUT_DIR/${s}.het.1kb" >/dev/null

  # append to merged table with the sample name
  awk -v S="$s" 'NR>1{print S"\t"$1"\t"$2"\t"$3"\t"$4}' \
      "$OUT_DIR/${s}.het.1kb.snpden" >> "$MERGED"
done

echo " Done. Results in: $OUT_DIR"
echo "   - Per sample:  ${s}.het.1kb.snpden  (CHROM  BIN_START  N_HET  SNP_DENSITY)"
echo "   - Merged:      $(basename "$MERGED")"

