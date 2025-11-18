
res <- run_coloc_locuscompare(
  gwas_file = "extract_results/around_rs12151021_ABCA7_gwas.tsv",
  eqtl_file = "extract_results/around_rs12151021_ABCA7_eqtl.tsv",
  fig_file  = "figures/locuscompare_rs12151021_ABCA7.png",
  title1    = "GWAS",
  title2    = "ACE eQTL (cortex-EUR)"
)

res$top_h4


# --------------- function for multiple candidates -------------- #

# library(data.table)
# library(coloc)
# library(locuscomparer)
# library(ggplot2)

run_coloc_locuscompare <- function(
    gwas_file,
    eqtl_file,
    fig_file,
    title1 = "GWAS",
    title2 = "eQTL",
    gwas_N_cases = 53042,      # 53,042 affected
    gwas_N_total = 408942,     # 53,042 cases + 355,900 controls
    sdY = 1,      # Spearman + z-meta
    out_prefix = NULL,    
    width = 16,
    height = 8,
    dpi = 300
) {
  # GWAS
  gwas <- fread(
    gwas_file,
    sep    = "\t",
    header = TRUE,
    select = c("hm_rsid", "hm_other_allele", "hm_effect_allele",
               "hm_beta", "standard_error",
               "hm_effect_allele_frequency", "p_value")
  )
  
  # eQTL
  eqtl <- fread(
    eqtl_file,
    sep    = "\t",
    header = TRUE,
    select = c("SNP", "SNPAlleles", "SNPEffectAllele",
               "SNPEffectAlleleFreq", "MetaP",
               "MetaPN", "MetaPZ", "MetaBeta", "MetaSE", "NrDatasets")
  )
  eqtl[, rsid := tstrsplit(SNP, ":", keep = 3)]
  
  # merge
  input <- merge(
    gwas,
    eqtl,
    by.x     = "hm_rsid",
    by.y     = "rsid",
    all      = FALSE,
    suffixes = c("_gwas", "_eqtl")
  )
  
  # flip
  input_flipped <- copy(input)
  
  idx_same <- input_flipped$SNPEffectAllele == input_flipped$hm_effect_allele
  idx_flip <- input_flipped$SNPEffectAllele == input_flipped$hm_other_allele
  idx_bad  <- !(idx_same | idx_flip)
  
  message(
    "Allele alignment: same = ", sum(idx_same),
    ", flip = ", sum(idx_flip),
    ", bad = ", sum(idx_bad)
  )
  
  bad_rows <- input[idx_bad]
  print(bad_rows)
  
  input_flipped[idx_flip, `:=`(
    SNPEffectAllele = hm_effect_allele,
    SNPEffectAlleleFreq = 1 - SNPEffectAlleleFreq,
    MetaBeta = -MetaBeta
  )]
  
  input_flipped <- input[!idx_bad]
  
  gwas_N <- gwas_N_cases
  gwas_s <- gwas_N_cases / gwas_N_total
  
  N_eqtl <- mean(input_flipped$MetaPN, na.rm = TRUE)
  
  result_flipped <- coloc.abf(
    dataset1 = list(
      snp = input_flipped$hm_rsid,
      pvalues = input_flipped$p_value,
      beta = input_flipped$hm_beta,
      varbeta = input_flipped$standard_error^2,
      type = "cc",
      N = gwas_N,
      s = gwas_s
    ),
    dataset2 = list(
      snp = input_flipped$hm_rsid,
      pvalues = input_flipped$MetaP,
      beta = input_flipped$MetaBeta,
      varbeta = input_flipped$MetaSE^2,
      type = "quant",
      N = N_eqtl,
      MAF = pmin(input_flipped$SNPEffectAlleleFreq, 1 - input_flipped$SNPEffectAlleleFreq),
      sdY = sdY
    )
  )
  
  res_tab <- result_flipped$results
  ord <- order(-res_tab$SNP.PP.H4)
  top_h4 <- head(res_tab[ord, ], 5)
  
  # locuscompare
  gwas_plot <- data.frame(
    rsid = gwas$hm_rsid,
    pval = gwas$p_value
  )
  eqtl_plot <- data.frame(
    rsid = eqtl$rsid,
    pval = eqtl$MetaP
  )
  
  gwas_plot$pval <- as.numeric(gwas_plot$pval)
  eqtl_plot$pval <- as.numeric(eqtl_plot$pval)
  
  p <- locuscompare(
    in_fn1 = gwas_plot,
    in_fn2 = eqtl_plot,
    title1 = title1,
    title2 = title2
  )
  
  ggsave(
    filename = fig_file,
    plot = p,
    width = width,
    height = height,
    dpi = dpi
  )

  if (is.null(out_prefix)) {
    base <- tools::file_path_sans_ext(basename(fig_file))
    out_prefix <- file.path(dirname(fig_file), base)
  }
  
  rds_file   <- paste0(out_prefix, "_coloc_result_flipped.rds")
  table_file <- paste0(out_prefix, "_coloc_snps_flipped.tsv")
  
  saveRDS(result_flipped, file = rds_file)
  fwrite(result_flipped$results, file = table_file, sep = "\t")
  
  message("Saved coloc result to: ", rds_file, " and ", table_file)
  message("Saved figure to: ", fig_file)
  
  invisible(list(
    result_flipped = result_flipped,
    top_h4 = top_h4,
    input_flipped  = input_flipped,
    plot = p
  ))
}
