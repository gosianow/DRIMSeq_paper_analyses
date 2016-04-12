######################################################
## <<geuvadis_drimseq_0_3_3_run_speed.R>>

# BioC 3.2
# Created 7 Feb 2016

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
# workers=5
# population='CEU'
# chr='19'

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

out_dir <- "drimseq_0_3_3_analysis_speed/"
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

d <- dmTest(d, BPPARAM = BPPARAM)


plotTest(d, out_dir = out_name)


save(d, file = paste0(out_name, "d.Rdata"))


res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

table(res$adj_pvalue < 0.05)


write.table(res, file = paste0(out_name, "results.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)







### Plot dispersion versus mean with marked significant sqtls
res_sign <- res[res$adj_pvalue < 0.05, , drop = FALSE]

ggp <- plotDispersion(d)

ggp2 <- ggp +
  geom_point(data = ggp$data[paste0(res_sign$gene_id, ".", res_sign$block_id), ], aes(x = mean_expression, y = dispersion), color = "black", size = 0.6)

pdf(paste0(out_name, "dispersion_vs_mean_marked_sqtls.pdf"))
print(ggp2)
dev.off()





### Plot dispersion versus mean with theta as dispersion

genewise_disp_theta <- 1 / (1 + 10^ggp$data$dispersion)

df <- data.frame(mean_expression = ggp$data$mean_expression, dispersion = genewise_disp_theta, nr_features = ggp$data$nr_features)

df_quant <- min(quantile(na.omit(df$nr_features), probs = 0.95), 30)
breaks <- seq(2, df_quant, ceiling(df_quant/10))


ggp2 <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion", colour = "nr_features" )) +
  theme_bw() +
  xlab("Log10 of mean expression") +
  ylab("Theta") +
  geom_point(alpha = 0.7, na.rm = TRUE) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") +
  guides(colour = guide_colorbar(barwidth = 20, barheight = 0.5)) +
  scale_colour_gradient(limits = c(2, max(breaks)), breaks = breaks, low = "royalblue2", high="red2", name = "Number of features", na.value = "red2")


pdf(paste0(out_name, "dispersion_vs_mean_theta.pdf"))
print(ggp2)
dev.off()





sessionInfo()












