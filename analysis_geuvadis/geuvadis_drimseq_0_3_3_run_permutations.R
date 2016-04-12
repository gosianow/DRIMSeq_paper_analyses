######################################################
## <<geuvadis_drimseq_0_3_3_run_permutations.R>>

# BioC 3.2 - R32dev
# Created 7 Feb 2016
# Modified 7 Apr 2016

##############################################################################
Sys.time()
##############################################################################

library(BiocParallel)
library(DRIMSeq)
library(ggplot2)
library(limma)
library(GenomicRanges)

##############################################################################
# Arguments for testing the code
##############################################################################

# rwd='/home/Shared/data/seq/geuvadis'
# workers=10
# population='CEU'
# chr='17'

##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

print(rwd)
print(workers)
print(chr)


##############################################################################

setwd(rwd)

out_dir <- "drimseq_0_3_3_analysis_permutations_all_genes/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_name <- paste0(out_dir, population, "_chr",chr, "_")

data_dir <- "data/"


if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}

########################################################
# sqtl analysis per chromosome
########################################################

### Input files: transcript expression, gene location and genotype information

### read genotypes
genotypes_raw <- read.table(paste0(data_dir, "genotypes/snps_", population, "_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(paste0("chr", chr), nrow(genotypes_raw)), IRanges::IRanges(genotypes_raw$start, genotypes_raw$end))
names(snp_ranges) <- genotypes_raw$snpId


### read ranges

genes_path <- paste0(data_dir, "annotation/gencode.v12.annotation.gtf")
gtf0 <- rtracklayer::import(genes_path)

## keep protein coding genes
keep <- mcols(gtf0)$gene_type == "protein_coding" & mcols(gtf0)$type == "gene" & seqnames(gtf0) == paste0("chr", chr)
gene_ranges <- gtf0[keep]

names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$gene_id


### read counts 
counts_path <- paste0(data_dir, "expression/trExpCount_", population, ".tsv")
counts_raw <- read.table(counts_path, header = TRUE, as.is = TRUE)

counts_raw <- counts_raw[counts_raw$geneId %in% names(gene_ranges), ]

stopifnot(all(strsplit2(colnames(counts_raw[, -c(1:2)]), "\\.")[, 1] == colnames(genotypes_raw[, -c(1:4)])))



### DRIMSeq SQTL analysis

d <- dmSQTLdataFromRanges(counts = counts_raw[, -c(1:2)], gene_id = counts_raw$geneId, feature_id = counts_raw$trId, gene_ranges = gene_ranges, genotypes = genotypes_raw[, -c(1:4)], snp_id = genotypes_raw$snpId, snp_ranges = snp_ranges, sample_id = colnames(genotypes_raw[, -c(1:4)]), window = 5e3, BPPARAM = BPPARAM)


rm("counts_raw", "genotypes_raw", "gene_ranges", "snp_ranges", "gtf0")


d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5, min_samps_feature_prop = 0, minor_allele_freq = 5, min_gene_expr = 10, min_feature_expr = 10, min_feature_prop = 0, max_features = Inf, BPPARAM = BPPARAM)


plotData(d, out_dir = out_name)


d <- dmDispersion(d, common_dispersion = FALSE, disp_init = 10, verbose = TRUE,  speed = TRUE, BPPARAM = BPPARAM)

plotDispersion(d, out_dir = out_name)


d <- dmFit(d, BPPARAM = BPPARAM)


d <- dmTest(d, permutations = "all_genes", verbose = 1, BPPARAM = BPPARAM)


plotTest(d, out_dir = out_name)


save(d, file = paste0(out_name, "d.Rdata"))


res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

table(res$adj_pvalue < 0.05)


write.table(res, file = paste0(out_name, "results.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#############################################################################

### Some investigation of the results

#############################################################################

# ### Plot dispersion versus mean with marked significant sqtls
# res_sign <- res[res$adj_pvalue < 0.05, , drop = FALSE]
# 
# ggp <- plotDispersion(d)
# 
# ggp2 <- ggp +
#   geom_point(data = ggp$data[paste0(res_sign$gene_id, ".", res_sign$block_id), ], aes(x = mean_expression, y = dispersion), color = "black", size = 0.6)
# 
# pdf(paste0(out_name, "dispersion_vs_mean_marked_sqtls.pdf"))
# print(ggp2)
# dev.off()
# 
# 
# 
# ### Read in the results from not permutated analysis
# res_lr <- read.table(paste0("drimseq_0_3_3_analysis/", population, "_chr",chr, "_results.txt"), header = TRUE, as.is = TRUE)
# 
# 
# res_lr <- unique(res_lr[, !grepl("snp_id", colnames(res_lr))])
# res <- unique(res[, !grepl("snp_id", colnames(res))])
# 
# 
# ### Add to the results the standard LR p-values 
# res$pvalue_lr <- pchisq(res$lr, df = res$df , lower.tail = FALSE)
# res$adj_pvalue_lr <- p.adjust(res$pvalue_lr, method="BH")
# 
# 
# table(res$adj_pvalue < 0.05)
# 
# table(res$adj_pvalue_lr < 0.05)
# 
# table(res_lr$adj_pvalue < 0.05)
# 
# 
# ### Plot permutation adjusted p-values versus LR standard p-values
# 
# ggp <- ggplot(res, aes(x = pvalue_lr, y = pvalue)) +
#   geom_point(alpha = 0.6) +
#   geom_abline(intercept = 0, slope = 1, color = "orange")
# 
# 
# pdf(paste0(out_name, "pval_perm_vs_pval_lr.pdf"))
# print(ggp)
# dev.off()
# 
# 
# 
# ggp <- ggplot(res, aes(x = -log10(pvalue_lr), y = -log10(pvalue))) +
#   geom_point(alpha = 0.6) +
#   geom_abline(intercept = 0, slope = 1, color = "orange")
# 
# 
# pdf(paste0(out_name, "pval_perm_vs_pval_lr_log.pdf"))
# print(ggp)
# dev.off()
# 
# 
# ### Plot p-values when speed is TRUE and FALSE
# ggdf <- merge(res, res_lr, by = c("gene_id", "block_id"), sort = FALSE)
# 
# ggp <- ggplot(ggdf, aes(x = -log10(pvalue_lr), y = -log10(pvalue.y))) +
#   geom_point(alpha = 0.6) +
#   geom_abline(intercept = 0, slope = 1, color = "orange") +
#   xlab("-log10(p-value) speed = TRUE") +
#   ylab("-log10(p-value) speed = FALSE")
# 
# 
# pdf(paste0(out_name, "pval_nospeed_vs_pval_speed_log.pdf"))
# print(ggp)
# dev.off()
# 
# 
# 
# ### Plot a histogram of standard LR p-values
# ggp <- DRIMSeq:::dm_plotPvalues(res$pvalue_lr)
# 
# pdf(paste0(out_name, "hist_pvalues_lr.pdf"))
# print(ggp)
# dev.off()
# 
# 
# 
# ### Read in sqtlseeker results
# 
# 
# res_seeker <- read.table(paste0("sqtlseeker_2_1_analysis/results/CEU_results_all.txt"), header = TRUE, as.is = TRUE)
# 
# res_seeker <- res_seeker[grep(paste0("snp_", chr, "_"), res_seeker$snpId), ]
# 
# 
# df <- data.frame(pvalues = res_seeker[, "pv"])
# 
# ggp <- ggplot(df, aes_string(x = "pvalues")) +
#   theme_bw() +
#   xlab("p-values") +
#   ylab("Frequency") +
#   geom_histogram(binwidth = 0.01, fill = "deeppink4") +
#   theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=16, face="bold")) +
#   coord_cartesian(xlim = c(-0.02, 1.02)) +
#   geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(nrow(df), " tests       ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 3, size = 6)
# 
# 
# 
# pdf(paste0(out_name, "hist_pvalues_sqtlseeker.pdf"))
# print(ggp)
# dev.off()
# 
# 
# 
# res <- results(d)
# 
# df <- data.frame(pvalues = res[, "pvalue"])
# 
# ggp <- ggplot(df, aes_string(x = "pvalues")) +
#   theme_bw() +
#   xlab("p-values") +
#   ylab("Frequency") +
#   geom_histogram(binwidth = 0.01, fill = "deeppink4") +
#   theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=16, face="bold")) +
#   coord_cartesian(xlim = c(-0.02, 1.02)) +
#   geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(nrow(df), " tests       ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 3, size = 6)
# 
# 
# 
# pdf(paste0(out_name, "hist_pvalues_snps.pdf"))
# print(ggp)
# dev.off()
# 
# 
# 
# 
# ### Plot p-values for DRIMSeq and sqtlseeker
# ggdf <- merge(res, res_seeker, by.x = c("gene_id", "snp_id"), by.y = c("geneId", "snpId"), sort = FALSE)
# 
# ggp <- ggplot(ggdf, aes(x = -log10(pvalue), y = -log10(pv))) +
#   geom_point(alpha = 0.6) +
#   geom_abline(intercept = 0, slope = 1, color = "orange") +
#   xlab("-log10(p-value) DRIMSeq") +
#   ylab("-log10(p-value) sQTLseekeR")
# 
# 
# pdf(paste0(out_name, "pval_drimseq_vs_sqtlseeker.pdf"))
# print(ggp)
# dev.off()
# 









sessionInfo()












