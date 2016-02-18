######################################################
## ----- geuvadis_drimseq_0_3_3_random
## <<geuvadis_drimseq_0_3_3_random.R>>

# BioC 3.2
# Created 9 Feb 2016

##############################################################################
# Do the sqtl analysis for SNP randomly assigned to genes. This is a validation method. We want to see how many significant sQTLs will be found. If many, then the method works badly.

Sys.time()

##############################################################################

library(DRIMSeq)
library(ggplot2)
library(limma)
library(GenomicRanges)

##############################################################################
# Arguments for testing the code
##############################################################################

# rwd='/home/Shared/data/seq/geuvadis'
# workers=4
# population='CEU'
# chr='22'

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

out_dir <- "drimseq_0_3_3_analysis_random/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_name <- paste0(out_dir, population, "_chr",chr, "_")

data_dir <- "data/"

########################################################
# sqtl analysis per chromosome for randommly assignes SNPs genes
########################################################


### read genotypes
genotypes_raw <- read.table(paste0(data_dir, "genotypes/snps_", population, "_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)


### read ranges

genes_path <- paste0(data_dir, "annotation/gencode.v12.annotation.gtf")
gtf0 <- rtracklayer::import(genes_path)

## keep protein coding genes
keep <- mcols(gtf0)$gene_type == "protein_coding" & mcols(gtf0)$type == "gene" & seqnames(gtf0) == paste0("chr", chr)
gene_ranges <- gtf0[keep]

names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$gene_id

end(gene_ranges) <- start(gene_ranges) + 1


### randomly assign genes to SNPs by setting SNP loci equal to the starting position of a gene 
set.seed(123)
genotypes_raw$start <- start(gene_ranges)[sample(x = length(gene_ranges), nrow(genotypes_raw), replace = TRUE)]
genotypes_raw$end <- genotypes_raw$start


# write randomly assigned SNP so they can be used for sqtlseeker
dir.create(paste0(data_dir, "genotypes/random"), recursive = TRUE, showWarnings = FALSE)

write.table(genotypes_raw, paste0(data_dir, "genotypes/random/snps_", population, "_chr", chr ,"_random.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(paste0("chr", chr), nrow(genotypes_raw)), IRanges::IRanges(genotypes_raw$start, genotypes_raw$end))
names(snp_ranges) <- genotypes_raw$snpId



### read counts 
counts_path <- paste0(data_dir, "expression/trExpCount_", population, ".tsv")
counts_raw <- read.table(counts_path, header = TRUE, as.is = TRUE)

counts_raw <- counts_raw[counts_raw$geneId %in% names(gene_ranges), ]

stopifnot(all(strsplit2(colnames(counts_raw[, -c(1:2)]), "\\.")[, 1] == colnames(genotypes_raw[, -c(1:4)])))


### DRIMSeq SQTL analysis
# use a window = 1

d <- dmSQTLdataFromRanges(counts = counts_raw[, -c(1:2)], gene_id = counts_raw$geneId, feature_id = counts_raw$trId, gene_ranges = gene_ranges, genotypes = genotypes_raw[, -c(1:4)], snp_id = genotypes_raw$snpId, snp_ranges = snp_ranges, sample_id = colnames(genotypes_raw[, -c(1:4)]), window = 1, BPPARAM = BiocParallel::MulticoreParam(workers = workers))


d <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5, min_samps_feature_prop = 0, minor_allele_freq = 5, min_gene_expr = 10, min_feature_expr = 10, min_feature_prop = 0, max_features = Inf, BPPARAM = BiocParallel::MulticoreParam(workers = workers))


plotData(d, out_dir = out_name)


d <- dmDispersion(d, verbose = TRUE,  speed = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotDispersion(d, out_dir = out_name)


d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = workers), verbose = TRUE)

d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = workers), verbose = TRUE)


plotTest(d, out_dir = out_name)


save(d, file = paste0(out_name, "d.Rdata"))


res <- results(d)
res <- res[order(res$pvalue, decreasing = FALSE), ]

table(res$adj_pvalue < 0.05)


write.table(res, file = paste0(out_name, "results.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


sessionInfo()












