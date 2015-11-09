# BioC 3.1

# Created 2 Nov 2015

##############################################################################################################

setwd("/home/Shared/data/seq/geuvadis/")


library(DRIMSeq)
library(ggplot2)

out_dir <- "drimseq_0_3_1_analysis/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


########################################################
# sqtl data for chr1
########################################################

## Input files: transcript expression, gene location and genotype information
data_dir <- "data/"
chr <- "1"


### read counts 
counts_path <- paste0(data_dir, "expression/trExpCount_CEU.tsv")
counts_raw <- read.table(counts_path, header = TRUE, as.is = TRUE)

gene_id <- counts_raw$geneId
feature_id <- counts_raw$trId
counts <- counts_raw[, -c(1:2)]


### read ranges

genes_path = paste0(data_dir, "annotation/gencode.v12.annotation_genes.bed")
gene_ranges = rtracklayer::import(genes_path)
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name


### read genotypes
genotypes_raw <- read.table(paste0(data_dir, "genotypes/snps_CEU_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(chr, nrow(genotypes_raw)), IRanges::IRanges(genotypes_raw$start, genotypes_raw$end))
names(snp_ranges) <- genotypes_raw$snpId
snp_id <- genotypes_raw$snpId

genotypes <- genotypes_raw[, -c(1:4)]


sample_id <- colnames(genotypes)

window <- 5e3


### DM SQTL analysis

d <- dmSQTLdataFromRanges(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotData(d, paste0(out_dir, "chr",chr, "_filt0_"))

save(d, file = paste0(out_dir, "chr",chr, "_d.Rdata"))


load(paste0(out_dir, "chr",chr, "_d.Rdata"))


d_filt1 <- dmFilter(d, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.05, max_features = Inf, minor_allele_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotData(d_filt1, paste0(out_dir, "chr",chr, "_filt1_"))


d_filt2 <- dmFilter(d, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 10, min_feature_prop = 0.05, max_features = Inf, minor_allele_freq = 10, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotData(d_filt2, paste0(out_dir, "chr",chr, "_filt2_"))



d_filt3 <- dmFilter(d, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 20, min_feature_prop = 0.01, max_features = Inf, minor_allele_freq = 20, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotData(d_filt3, paste0(out_dir, "chr",chr, "_filt3_"))



d_filt4 <- dmFilter(d, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 20, min_feature_prop = 0.05, max_features = Inf, minor_allele_freq = 20, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotData(d_filt3, paste0(out_dir, "chr",chr, "_filt4_"))


# library(devtools)
# load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")
# 
# x = d
# mean_expression = TRUE
# common_dispersion = TRUE
# genewise_dispersion = TRUE
# disp_adjust = TRUE
# disp_mode = "grid"
# disp_interval = c(0, 1e+4)
# disp_tol = 1e-08
# disp_init = 100
# disp_init_weirMoM = TRUE
# disp_grid_length = 21
# disp_grid_range = c(-10, 10)
# disp_moderation = "none"
# disp_prior_df = 1
# disp_span = 0.3
# prop_mode = "constrOptimG"
# prop_tol = 1e-12
# verbose = TRUE
# speed = TRUE
# BPPARAM = BiocParallel::MulticoreParam(workers = 10)
# 
# counts = x@counts
# disp_tol = 1e+01
# 
# gamma0 = 3819.66
# 
# dispersion = gamma0
# model = "full"



d <- dmDispersion(d_filt4, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 10), speed = FALSE)

plotDispersion(d, paste0(out_dir, "chr",chr, "_filt4b_"))


d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = 10))

plotTest(d, out_dir = paste0(out_dir, "chr",chr, "_filt4b_"))


save(d, file = paste0(out_dir, "chr",chr, "_filt4b_d_dmSQTLtest.Rdata"))

res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

table(res$adj_pvalue < 0.05)

gene_id <- res$gene_id[1:10]
snp_id <- res$snp_id[1:10]


plotFit(d, gene_id, snp_id, out_dir = paste0(out_dir, "chr",chr, "_filt4b_"), plot_type = "boxplot1")























