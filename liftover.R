
# install packages

install.packages("data.table")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("rtracklayer", "GenomicRanges"))

install.packages('R.utils')

# load packages

library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(R.utils)

# read GWAS harmonised file

gwas <- fread("33589840-GCST90012877-EFO_0000249.h.tsv.gz")

# check columns
colnames(gwas)

# read UCSC chain file
# chain <- import.chain("hg19ToHg38.over.chain.gz") Error why?
chain <- import.chain("hg19ToHg38.over.chain")

# check if NAs in the GWAS harmonised file
sum(is.na(gwas$hm_pos)) # 41093
sum(is.na(gwas$hm_chrom)) # 41093

# remove NAs in the GWAS harmonised file
gwas_noNA <- gwas[!is.na(hm_chrom) & !is.na(hm_pos)]

nrow(gwas) # 10648365
nrow(gwas_noNA) # 10607272

# convert hg19 to GRanges
gr_hg19 <- GRanges(
  seqnames = paste0("chr", gwas_noNA$hm_chrom),
  ranges   = IRanges(start = gwas_noNA$hm_pos, width = 1),
  hm_rsid  = gwas_noNA$hm_rsid
)

# liftover: hg19 to hg38
gr_hg38_list <- liftOver(gr_hg19, chain)

# only remain SNP width = 1bp
idx_1to1 <- elementNROWS(gr_hg38_list) == 1
gr_hg38_1to1 <- gr_hg38_list[idx_1to1]

# generate GRanges
gr_hg38 <- unlist(gr_hg38_1to1)

# filter again: width = 1bp
gr_hg38 <- gr_hg38[width(gr_hg38) == 1]

# from GRanges to hg38
map_hg38 <- data.table(
  hm_rsid  = mcols(gr_hg38)$hm_rsid,
  chr_hg38 = as.character(seqnames(gr_hg38)),
  pos_hg38 = start(gr_hg38)  # new coordinate
)

# remove "chr"
map_hg38[, chr_hg38 := sub("^chr", "", chr_hg38)]

# merge hg38 and GWAS
gwas_hg38 <- merge(
  gwas,
  map_hg38,
  by = "hm_rsid",
  all.x = TRUE
)

fwrite(gwas_hg38,
       file = "GWAS_hg38.tsv.gz",
       sep = "\t")


