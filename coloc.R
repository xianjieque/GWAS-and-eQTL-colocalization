
install.packages("coloc")
install.packages("devtools")
devtools::install_github("boxiangliu/locuscomparer")

library(coloc)
library(locuscomparer)

gwas <- fread("extract_results/around_rs4311_ACE_gwas.tsv", 
              sep = "\t", 
              header = T, 
              select = c("hm_rsid", "hm_other_allele", "hm_effect_allele", 
                         "hm_beta", "standard_error", "hm_effect_allele_frequency", "p_value"))
head(gwas)

eqtl <- fread("extract_results/around_rs4311_ACE_eqtl.tsv", 
              sep = "\t", 
              header = T, 
              select = c("SNP", "SNPAlleles", "SNPEffectAllele", 
                         "SNPEffectAlleleFreq", "MetaP", 
                         "MetaPN", "MetaPZ", "MetaBeta", "MetaSE", "NrDatasets"))
eqtl$rsid <- sapply(strsplit(eqtl$SNP, ":"), `[`, 3)
head(eqtl)

# merge, find overlap SNP
input <- merge(gwas, eqtl, 
               by.x = "hm_rsid", 
               by.y = "rsid", 
               all = F, 
               suffixes = c("_gwas", "_eqtl"))
head(input)

# colocalize
result <- coloc.abf(
  dataset1 = list(snp = input$hm_rsid, 
                  pvalues = input$p_value, 
                  beta = input$hm_beta, 
                  varbeta = input$standard_error^2, 
                  type = "cc",  # case-control GWAS
                  N = 53042,
                  s = 53042/408942), 
  dataset2 = list(snp = input$hm_rsid, 
                  pvalues = input$MetaP, 
                  beta = input$MetaBeta, 
                  varbeta = input$MetaSE^2, 
                  type = "quant", 
                  N = mean(eqtl$MetaPN), 
                  MAF = pmin(input$SNPEffectAlleleFreq, 1 - input$SNPEffectAlleleFreq),
                  sdY = 1),
)

subset(result$results, SNP.PP.H4 > 0.5)

# here I set "N = 53042, s = 53042/408942" because "There were 898 AD cases, 52,791 AD proxy cases 
# and 355,900 controls in the combined white British and white non-British cohorts. 
# For the association analyses, we lumped the true and proxy cases together 
# (53,042 unique affected individuals) and used the linear mixed model implemented in BOLT-LMM93 v.2.3.2."

# I set "sdY = 1" because "We performed cis-eQTL analysis in each of the eQTL discovery datasets by 
# calculating Spearman correlations within each cohort, followed by a sample size-weighted z-score meta-analysis approach"

# In the merged dataset "input", in some rows, GWAS's hm_effect_allele and eQTL's effect allele "SNPEffectAllele"
# are exactly reversed. It won't influence the colocalization results but will influence interpretation.

# ---------------flip-------------- #

input_flipped <- copy(input) 

idx_same  <- input_flipped$SNPEffectAllele == input_flipped$hm_effect_allele

idx_flip  <- input_flipped$SNPEffectAllele == input_flipped$hm_other_allele

idx_bad   <- !(idx_same | idx_flip)

cat("same:", sum(idx_same),
    " flip:", sum(idx_flip),
    " bad:", sum(idx_bad), "\n")

# same: 1538  flip: 211  bad: 0 

input_flipped[idx_flip, `:=`(
  SNPEffectAllele = hm_effect_allele,
  SNPEffectAlleleFreq = 1 - SNPEffectAlleleFreq,
  MetaBeta = -MetaBeta
)]

# input_flipped <- input_flipped[!idx_bad]

result_flipped <- coloc.abf(
  dataset1 = list(snp = input_flipped$hm_rsid, 
                  pvalues = input_flipped$p_value, 
                  beta = input_flipped$hm_beta, 
                  varbeta = input_flipped$standard_error^2, 
                  type = "cc",  # case-control GWAS
                  N = 53042,
                  s = 53042/408942), 
  dataset2 = list(snp = input_flipped$hm_rsid, 
                  pvalues = input_flipped$MetaP, 
                  beta = input_flipped$MetaBeta, 
                  varbeta = input_flipped$MetaSE^2, 
                  type = "quant", 
                  N = mean(eqtl$MetaPN), 
                  MAF = pmin(input_flipped$SNPEffectAlleleFreq, 1 - input_flipped$SNPEffectAlleleFreq),
                  sdY = 1),
)

subset(result_flipped$results, SNP.PP.H4 > 0.5)

# --------------- visualize -------------- #

fwrite(eqtl, "extract_results/around_rs4311_ACE_eqtl_with_rsid.tsv", sep = "\t")

p <- locuscompare(
  in_fn1 = "extract_results/around_rs4311_ACE_gwas.tsv", 
  in_fn2 = "extract_results/around_rs4311_ACE_eqtl_with_rsid.tsv", 
  title1 = "GWAS", 
  title2 = "ACE eQTL (cortex-EUR)", 
  marker_col1 = "hm_rsid", 
  pval_col1 = "p_value", 
  marker_col2 = "rsid", 
  pval_col2   = "MetaP"
  )

ggsave("figures/locuscompare_rs4311_ACE.png", p, width=8, height=8, dpi=300)

