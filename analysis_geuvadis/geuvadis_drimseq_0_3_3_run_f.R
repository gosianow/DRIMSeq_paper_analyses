######################################################
## ----- geuvadis_drimseq_0_3_3_run_f
## <<geuvadis_drimseq_0_3_3_run_f.R>>

# BioC 3.2 - R32dev
# Created 10 Feb 2016

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

out_dir <- "drimseq_0_3_3_analysis_f/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_name <- paste0(out_dir, population, "_chr",chr, "_")

data_dir <- "data/"

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

d <- dmSQTLdataFromRanges(counts = counts_raw[, -c(1:2)], gene_id = counts_raw$geneId, feature_id = counts_raw$trId, gene_ranges = gene_ranges, genotypes = genotypes_raw[, -c(1:4)], snp_id = genotypes_raw$snpId, snp_ranges = snp_ranges, sample_id = colnames(genotypes_raw[, -c(1:4)]), window = 5e3, BPPARAM = BiocParallel::MulticoreParam(workers = workers))

rm("counts_raw", "genotypes_raw", "gene_ranges", "snp_ranges")


do <- dmFilter(d, min_samps_gene_expr = 70, min_samps_feature_expr = 5, min_samps_feature_prop = 0, minor_allele_freq = 5, min_gene_expr = 10, min_feature_expr = 10, min_feature_prop = 0, max_features = Inf, BPPARAM = BiocParallel::MulticoreParam(workers = workers))


plotData(do, out_dir = out_name)

rm("d")



################################################
### With CR adjustment
################################################

# d <- dmDispersion(do, verbose = TRUE,  speed = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# 
# plotDispersion(d, out_dir = paste0(out_name, "cr_"))
# 
# 
# d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# 
# ## LR test
# d <- dmTest(d, test = "lr", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# plotTest(d, out_dir = paste0(out_name, "lr_cr_"))
# 
# res <- results(d)
# 
# write.table(res, file = paste0(out_name, "results_lr_cr.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# 
# res <- res[order(res$pvalue, decreasing = FALSE), ]
# 
# plotFit(d, gene_id = res[1, "gene_id"], snp_id = res[1, "snp_id"], out_dir = paste0(out_name, "lr_cr_"))
# 
# 
# save(d, file = paste0(out_name, "d_cr.Rdata"))
# 
# ## F test
# d <- dmTest(d, test = "f", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# plotTest(d, out_dir = paste0(out_name, "f_cr_"))
# 
# res <- results(d)
# 
# write.table(res, file = paste0(out_name, "results_f_cr.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# 
# res <- res[order(res$pvalue, decreasing = FALSE), ]
# 
# plotFit(d, gene_id = res[1, "gene_id"], snp_id = res[1, "snp_id"], out_dir = paste0(out_name, "f_cr_"))
# 
# 

# 
# rm("d")
# 
# 
# 
# ################################################
# ### With no CR adjustment
# ################################################
# 
# 
# d <- dmDispersion(do, disp_adjust = FALSE, verbose = TRUE,  speed = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# 
# plotDispersion(d, out_dir = paste0(out_name, "ncr_"))
# 
# 
# d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# 
# ## LR test
# d <- dmTest(d, test = "lr", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# plotTest(d, out_dir = paste0(out_name, "lr_ncr_"))
# 
# res <- results(d)
# 
# write.table(res, file = paste0(out_name, "results_lr_ncr.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# 
# res <- res[order(res$pvalue, decreasing = FALSE), ]
# 
# plotFit(d, gene_id = res[1, "gene_id"], snp_id = res[1, "snp_id"], out_dir = paste0(out_name, "lr_ncr_"))
# 
# 
# save(d, file = paste0(out_name, "d_ncr.Rdata"))
# 
# ## F test
# d <- dmTest(d, test = "f", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
# 
# plotTest(d, out_dir = paste0(out_name, "f_ncr_"))
# 
# res <- results(d)
# 
# write.table(res, file = paste0(out_name, "results_f_ncr.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
# 
# 
# res <- res[order(res$pvalue, decreasing = FALSE), ]
# 
# plotFit(d, gene_id = res[1, "gene_id"], snp_id = res[1, "snp_id"], out_dir = paste0(out_name, "f_ncr_"))
# 
# 

# 
# rm("d")

################################################
### With CR adjustment and moderation
################################################

d <- dmDispersion(do, disp_adjust = TRUE, disp_moderation = "common", verbose = TRUE,  speed = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))


plotDispersion(d, out_dir = paste0(out_name, "crm_"))


d <- dmFit(d, BPPARAM = BiocParallel::MulticoreParam(workers = workers))


## LR test
d <- dmTest(d, test = "lr", BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotTest(d, out_dir = paste0(out_name, "lr_crm_"))

res <- results(d)

write.table(res, file = paste0(out_name, "results_lr_crm.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


res <- res[order(res$pvalue, decreasing = FALSE), ]

plotFit(d, gene_id = res[1, "gene_id"], snp_id = res[1, "snp_id"], out_dir = paste0(out_name, "lr_crm_"))


save(d, file = paste0(out_name, "d_crm.Rdata"))


## F test
d <- dmTest(d, test = "f", BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotTest(d, out_dir = paste0(out_name, "f_crm_"))

res <- results(d)

write.table(res, file = paste0(out_name, "results_f_crm.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


res <- res[order(res$pvalue, decreasing = FALSE), ]

plotFit(d, gene_id = res[1, "gene_id"], snp_id = res[1, "snp_id"], out_dir = paste0(out_name, "f_crm_"))


rm("d")



sessionInfo()












