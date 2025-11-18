#!/bin/bash

GWAS_FILE="GWAS_stats_hg19/33589840-GCST90012877-EFO_0000249.h.tsv.gz"
EQTL_FILE="eQTL_stats_hg38/2021-07-23-cortex-EUR-80PCs-chr17.txt.gz"
Variant="rs4311"
Gene="ACE"
OUTFILE="extract_results/around_${Variant}_gwas.tsv"
OUTFILE2="extract_results/around_${Variant}_eqtl.tsv"

# Find the chromosome and coordiante of the SNP
read CHR POS < <(zcat "${GWAS_FILE}" | awk -F'\t' -v lead="${Variant}" '$2 == lead {print $3, $4; exit}')
echo "FOUND AT ${CHR}:${POS}"

# window = e6
START=$((POS - 500000))
END=$((POS + 500000))

# Subset the SNPs within 50,000bp
{
  zcat "${GWAS_FILE}" | head -n 1
  zcat "${GWAS_FILE}" | \
    awk -F'\t' -v chr="${CHR}" -v s="${START}" -v e="${END}" '$3 == chr && $4 >= s && $4 <= e'
} > "${OUTFILE}"

echo "Processing eQTL"

# Subset the SNPs within 50,000bp
{
  zcat "${EQTL_FILE}" | head -n 1
  zcat "${EQTL_FILE}" | \
    awk -F'\t' -v chr="${CHR}" -v s="${START}" -v e="${END}" -v gene="${Gene}" '$5 == gene && $7 == chr && $8 >= s && $8 <= e {print}'
} > ${OUTFILE2}
