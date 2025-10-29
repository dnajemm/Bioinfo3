#!/usr/bin/env bash
set -euo pipefail

IN_DIR="/data/projet3/Output/results_1kb"          # contient *.het.1kb.snpden
OUT_DIR="/data/projet3/Output/results_windows"     # sorties 5/10/20 kb
mkdir -p "$OUT_DIR"

# --- fonction AWK : agrège N bins de 1 kb en fenêtres non-chevauchantes ---
aggregate_N() {
  local infile="$1"   # ex: YJS8112.het.1kb.snpden
  local N="$2"        # 5 -> 5 kb, 10 -> 10 kb, 20 -> 20 kb
  local outfile="$3"  # ex: YJS8112.het.5kb.tsv

  awk -v N="$N" 'BEGIN{OFS="\t"}
  NR==1 {next}                              # saute l’entête vcftools
  {
    c=$1; s=$2; n=$3;                      # CHROM BIN_START N_SNPS [...]
    if (!init){cur=c; init=1}
    if (c!=cur){                           # changement de chromosome -> reset
      i=0; cur=c
    }
    i++; pos[i]=s; val[i]=n
    if (i==N){                             # quand on a N bins de 1 kb
      start=pos[1]; end=s+1000             # fenêtre exacte de N*1kb
      sum=0; for (j=1;j<=N;j++) sum+=val[j]
      print c, start, end, sum             # CHROM START END N_HET
      i=0                                  # reset pour la fenêtre suivante
    }
  }' "$infile" > "$outfile"
}

# fichiers fusionnés (un par taille)
MERGED5="$OUT_DIR/all_samples.het.5kb.tsv"
MERGED10="$OUT_DIR/all_samples.het.10kb.tsv"
MERGED20="$OUT_DIR/all_samples.het.20kb.tsv"
echo -e "SAMPLE\tCHROM\tSTART\tEND\tN_HET\tWINDOW_BP"  > "$MERGED5"
echo -e "SAMPLE\tCHROM\tSTART\tEND\tN_HET\tWINDOW_BP"  > "$MERGED10"
echo -e "SAMPLE\tCHROM\tSTART\tEND\tN_HET\tWINDOW_BP"  > "$MERGED20"

for f in "$IN_DIR"/*.het.1kb.snpden; do
  s=$(basename "$f" .het.1kb.snpden)                 # ex: YJS8112

  # 5 kb
  aggregate_N "$f" 5  "$OUT_DIR/${s}.het.5kb.tsv"
  awk -v S="$s" -v W=5000  'NR>0{print S,$1,$2,$3,$4,W}' OFS="\t" \
      "$OUT_DIR/${s}.het.5kb.tsv"  >> "$MERGED5"

  # 10 kb
  aggregate_N "$f" 10 "$OUT_DIR/${s}.het.10kb.tsv"
  awk -v S="$s" -v W=10000 'NR>0{print S,$1,$2,$3,$4,W}' OFS="\t" \
      "$OUT_DIR/${s}.het.10kb.tsv" >> "$MERGED10"

  # 20 kb
  aggregate_N "$f" 20 "$OUT_DIR/${s}.het.20kb.tsv"
  awk -v S="$s" -v W=20000 'NR>0{print S,$1,$2,$3,$4,W}' OFS="\t" \
      "$OUT_DIR/${s}.het.20kb.tsv" >> "$MERGED20"
done

echo "OK. Sorties dans $OUT_DIR"
echo "   - Par échantillon : YJSxxxx.het.{5kb,10kb,20kb}.tsv  (CHROM START END N_HET)"
echo "   - Fusionnés       : $(basename "$MERGED5"), $(basename "$MERGED10"), $(basename "$MERGED20")"
