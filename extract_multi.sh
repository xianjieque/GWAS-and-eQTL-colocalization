#!/bin/bash

GWAS_FILE="GWAS_stats_hg19/33589840-GCST90012877-EFO_0000249.h.tsv.gz"
TARGET_LIST="targets.txt"   # row names: rsID GeneSymbol
WINDOW=500000               # 1 Megabase

mkdir -p extract_results

while read -r Variant Gene || [ -n "$Variant" ]; do
  [[ -z "$Variant" ]] && continue
  [[ "$Variant" =~ ^# ]] && continue

  echo "Processing Variant=${Variant}, Gene=${Gene}"

  read CHR POS < <(
    zcat "${GWAS_FILE}" \
      | awk -F'\t' -v lead="${Variant}" '$2 == lead {print $3, $4; exit}'
  )

  echo "  FOUND GWAS lead at ${CHR}:${POS}"

  START=$((POS - WINDOW))
  END=$((POS + WINDOW))

  EQTL_FILE="eQTL_stats_hg38/2021-07-23-cortex-EUR-80PCs-chr${CHR}.txt.gz"

  OUT_GWAS="extract_results/around_${Variant}_${Gene}_gwas.tsv"
  OUT_EQTL="extract_results/around_${Variant}_${Gene}_eqtl.tsv"

  {
    zcat "${GWAS_FILE}" | head -n 1
    zcat "${GWAS_FILE}" \
      | awk -F'\t' -v chr="${CHR}" -v s="${START}" -v e="${END}" '
          $3 == chr && $4 >= s && $4 <= e
        '
  } > "${OUT_GWAS}"

  echo "  GWAS window written to ${OUT_GWAS}"

  {
    zcat "${EQTL_FILE}" | head -n 1
    zcat "${EQTL_FILE}" \
      | awk -F'\t' -v chr="${CHR}" -v s="${START}" -v e="${END}" -v gene="${Gene}" '
          $5 == gene && $7 == chr && $8 >= s && $8 <= e
        '
  } > "${OUT_EQTL}"

done < "${TARGET_LIST}"
