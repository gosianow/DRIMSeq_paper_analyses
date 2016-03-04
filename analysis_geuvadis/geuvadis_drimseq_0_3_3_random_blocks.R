######################################################
## ----- geuvadis_drimseq_0_3_3_random_blocks
## <<geuvadis_drimseq_0_3_3_random_blocks.R>>

# BioC 3.2
# Created 29 Feb 2016

##############################################################################
# Do the sqtl analysis for BLOCKS (unique SNPs) randomly assigned to genes. This is a validation method. We want to see how many significant sQTLs will be found. If many, then the method works badly.

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

out_dir <- "drimseq_0_3_3_analysis_random_blocks/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_name <- paste0(out_dir, population, "_chr",chr, "_")

data_dir <- "data/"

########################################################
# sqtl analysis per chromosome for randommly assignes BLOCKS genes
########################################################


### read genotypes
genotypes_raw <- read.table(paste0(data_dir, "genotypes/snps_", population, "_chr", chr ,".tsv"), header = TRUE, sep = "\t", as.is = TRUE)

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(paste0("chr", chr), nrow(genotypes_raw)), IRanges::IRanges(genotypes_raw$start, genotypes_raw$end))
names(snp_ranges) <- genotypes_raw$snpId



### read gene ranges

genes_path <- paste0(data_dir, "annotation/gencode.v12.annotation.gtf")
gtf0 <- rtracklayer::import(genes_path)

### keep protein coding genes
keep <- mcols(gtf0)$gene_type == "protein_coding" & mcols(gtf0)$type == "gene" & seqnames(gtf0) == paste0("chr", chr)
gene_ranges <- gtf0[keep]
names(gene_ranges) <- S4Vectors::mcols(gene_ranges)$gene_id



### read counts 
counts_path <- paste0(data_dir, "expression/trExpCount_", population, ".tsv")
counts_raw <- read.table(counts_path, header = TRUE, as.is = TRUE)

counts_raw <- counts_raw[counts_raw$geneId %in% names(gene_ranges), ]

stopifnot(all(strsplit2(colnames(counts_raw[, -c(1:2)]), "\\.")[, 1] == colnames(genotypes_raw[, -c(1:4)])))


### Find blocks of SNPs

d <- dmSQTLdataFromRanges(counts = counts_raw[, -c(1:2)], gene_id = counts_raw$geneId, feature_id = counts_raw$trId, gene_ranges = gene_ranges, genotypes = genotypes_raw[, -c(1:4)], snp_id = genotypes_raw$snpId, snp_ranges = snp_ranges, sample_id = colnames(genotypes_raw[, -c(1:4)]), window = 5e3, BPPARAM = BiocParallel::MulticoreParam(workers = workers))


genotypes <- d@genotypes@unlistData


blocks2snps <- unlist(lapply(1:length(d@blocks), function(i){
  # i = 1
  blocks <- d@blocks[[i]]
  snp_names <- blocks[!duplicated(blocks[, "block_id"]), "snp_id"]
  
  return(snp_names)
  }), use.names = FALSE)


rownames(genotypes) <- blocks2snps



### randomly assign genes to SNPs by setting SNP loci equal to the starting position of a gene 
set.seed(123)
snp_start <- start(gene_ranges)[sample(x = length(gene_ranges), nrow(genotypes), replace = TRUE)]

snp_ranges <- GenomicRanges::GRanges(S4Vectors::Rle(paste0("chr", chr), nrow(genotypes)), IRanges::IRanges(snp_start, snp_start))
names(snp_ranges) <- rownames(genotypes)


genotypes_new <- data.frame(chr = chr, start = start(snp_ranges), end = end(snp_ranges), snpId = names(snp_ranges), genotypes, row.names = NULL, stringsAsFactors = FALSE)


# write randomly assigned SNP so they can be used for sqtlseeker
dir.create(paste0(data_dir, "genotypes/random_blocks"), recursive = TRUE, showWarnings = FALSE)

write.table(genotypes_new, paste0(data_dir, "genotypes/random_blocks/snps_", population, "_chr", chr ,"_random_blocks.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



# Like this SNPs are uniquely assigned to genes
end(gene_ranges) <- start(gene_ranges) + 1





### DRIMSeq SQTL analysis
# use a window = 1

d <- dmSQTLdataFromRanges(counts = counts_raw[, -c(1:2)], gene_id = counts_raw$geneId, feature_id = counts_raw$trId, gene_ranges = gene_ranges, genotypes = genotypes_new[, -c(1:4)], snp_id = genotypes_new$snpId, snp_ranges = snp_ranges, sample_id = colnames(genotypes_new[, -c(1:4)]), window = 1, BPPARAM = BiocParallel::MulticoreParam(workers = workers))


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












