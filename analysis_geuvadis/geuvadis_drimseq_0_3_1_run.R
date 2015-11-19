######################################################
## ----- geuvadis_drimseq_0_3_1_run
## <<geuvadis_drimseq_0_3_1_run.R>>

# BioC 3.1
# Created 13 Nov 2015 

##############################################################################

library(DRIMSeq)
library(ggplot2)
library(limma)

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/Shared/data/seq/geuvadis'
workers=5
population='CEU'
chr='19'

##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(workers)
print(chr)


##############################################################################

setwd(rwd)

out_dir <- "drimseq_0_3_1_analysis/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_name <- paste0(out_dir, population, "_chr",chr, "_")

data_dir <- "data/"

########################################################
# sqtl analysis per chromosome
########################################################

### Input files: transcript expression, gene location and genotype information

### read counts 
counts_path <- paste0(data_dir, "expression/trExpCount_", population, ".tsv")
counts_raw <- read.table(counts_path, header = TRUE, as.is = TRUE)

gene_id <- counts_raw$geneId
feature_id <- counts_raw$trId
counts <- counts_raw[, -c(1:2)]


### read ranges

genes_path = paste0(data_dir, "annotation/gencode.v12.annotation_genes.bed")
gene_ranges = rtracklayer::import(genes_path)
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$name


### read genotypes
genotypes_raw <- read.table(paste0(data_dir, "genotypes/snps_", population, "_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(chr, nrow(genotypes_raw)), IRanges::IRanges(genotypes_raw$start, genotypes_raw$end))
names(snp_ranges) <- genotypes_raw$snpId
snp_id <- genotypes_raw$snpId

genotypes <- genotypes_raw[, -c(1:4)]


all(strsplit2(colnames(counts), "\\.")[, 1] == colnames(genotypes))

sample_id <- colnames(genotypes)

window <- 5e3


### DRIMSeq SQTL analysis

d <- dmSQTLdataFromRanges(counts, gene_id, feature_id, gene_ranges, genotypes, snp_id, snp_ranges, sample_id, window, BPPARAM = BiocParallel::MulticoreParam(workers = workers))

# save(d, file = paste0(out_name, "d.Rdata"))
# load(paste0(out_name, "d.Rdata"))


d <- dmFilter(d, min_samps_gene_expr = 70, min_gene_expr = 1, min_samps_feature_prop = 5, min_feature_prop = 0.05, max_features = Inf, minor_allele_freq = 5, BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotData(d, out_dir = out_name)


d <- dmDispersion(d, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = workers), speed = FALSE)

plotDispersion(d, out_dir = out_name)


d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = workers))

d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = workers))


plotTest(d, out_dir = out_name)


save(d, file = paste0(out_name, "d.Rdata"))


res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

table(res$adj_pvalue < 0.05)


write.table(res, file = paste0(out_name, "results.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


sessionInfo()












