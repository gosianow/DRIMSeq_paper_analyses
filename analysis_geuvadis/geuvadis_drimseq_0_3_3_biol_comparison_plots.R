##############################################################################
## <<geuvadis_drimseq_0_3_3_biol_comparison_plots.R>>

# BioC 3.2
# Created 21 Apr 2016
# Modified


##############################################################################
Sys.time()
##############################################################################
# Libraries
##############################################################################

library(ggplot2)
library(reshape2)
library(plyr)
library(rtracklayer)
library(GenomicRanges)
library(limma)
library(BiocParallel)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/geuvadis'
# population='CEU'
# comparison_out='drimseq_0_3_3_comparison_permutations_all_genes'
# FDR=0.05
# path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf'
# workers=10

##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)


##############################################################################

setwd(rwd)

comparison_out <- paste0(comparison_out, "/")
dir.create(comparison_out, showWarnings = FALSE, recursive = TRUE)

out_dir <- comparison_out


data_dir <- "data/"


if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}


##############################################################################
# colors
##############################################################################


load(paste0(rwd, "/", comparison_out, "/colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)

#######################################################
# merge results 
#######################################################

results <- list()
metadata <- list()

#####################################
### sqtlseeker results
#####################################


res <- read.table(paste0(comparison_out, "results_sqtlseeker.txt"), header = TRUE, as.is = TRUE)
head(res)


results[["sqtlseeker"]] <- res
metadata[["sqtlseeker"]] <- data.frame(method_name = "sqtlseeker", stringsAsFactors = FALSE)


#####################################
### DRIMSeq results 
#####################################

res <- read.table(paste0(comparison_out, "results_drimseq.txt"), header = TRUE, as.is = TRUE)
head(res)

results[["drimseq"]] <- res
metadata[["drimseq"]] <- data.frame(method_name = "drimseq", stringsAsFactors = FALSE)


#####################################
### merge
#####################################

results_padj <- list()

padj_tmp <- results[["sqtlseeker"]][, c("gene_snp", "adj_pvalue")]
colnames(padj_tmp) <- c("gene_snp", "sqtlseeker")
results_padj[["sqtlseeker"]] <- padj_tmp


padj_tmp <- results[["drimseq"]][, c("gene_snp", "adj_pvalue")]
colnames(padj_tmp) <- c("gene_snp", "drimseq")
results_padj[["drimseq"]] <- padj_tmp


results_padj <- Reduce(function(...) merge(..., by = "gene_snp", all=TRUE, sort = FALSE), results_padj)
rownames(results_padj) <- results_padj$gene_snp




#####################################
### Get the gene lists
#####################################


res_sign_sqtlseeker <- results[["sqtlseeker"]][results[["sqtlseeker"]]$adj_pvalue < FDR, , drop = FALSE]
res_sign_drimseq <- results[["drimseq"]][results[["drimseq"]]$adj_pvalue < FDR, , drop = FALSE]


all_genes_sqtlseeker <- unique(results[["sqtlseeker"]]$gene_id)
all_genes_drimseq <- unique(results[["drimseq"]]$gene_id)


all_genes <- union(results[["sqtlseeker"]]$gene_id, results[["drimseq"]]$gene_id)

genes_sign_sqtlseeker <- unique(res_sign_sqtlseeker$gene_id)
genes_sign_drimseq <- unique(res_sign_drimseq$gene_id)

genes_sign_overlap <- intersect(genes_sign_sqtlseeker, genes_sign_drimseq)
genes_sign_sqtlseeker_unique <- setdiff(genes_sign_sqtlseeker, genes_sign_overlap)
genes_sign_drimseq_unique <- setdiff(genes_sign_drimseq, genes_sign_overlap)



all_sqtls <- union(results[["sqtlseeker"]]$gene_snp, results[["drimseq"]]$gene_snp)

sqtls_sign_sqtlseeker <- unique(res_sign_sqtlseeker$gene_snp)
sqtls_sign_drimseq <- unique(res_sign_drimseq$gene_snp)

sqtls_sign_overlap <- intersect(sqtls_sign_sqtlseeker, sqtls_sign_drimseq)
sqtls_sign_sqtlseeker_unique <- setdiff(sqtls_sign_sqtlseeker, sqtls_sign_overlap)
sqtls_sign_drimseq_unique <- setdiff(sqtls_sign_drimseq, sqtls_sign_overlap)



############################################################################
# Check the mean gene expression for the discovered sQTLs
# Check the number of expressed transcripts
############################################################################


### read counts 
counts_path <- paste0(data_dir, "expression/trExpCount_", population, ".tsv")
counts_raw <- read.table(counts_path, header = TRUE, as.is = TRUE)

counts <- counts_raw[, -grep("trId|geneId", colnames(counts_raw))]



### calculate mean gene expression
mean_expression <- by(counts, factor(counts_raw$geneId), function(x){
  mean(colSums(x, na.rm = TRUE), na.rm = TRUE)
}, simplify = FALSE)

mean_expression <- unlist(mean_expression)


ggdf <- data.frame(mean_expression = c(mean_expression[all_genes], mean_expression[genes_sign_sqtlseeker_unique], mean_expression[genes_sign_drimseq_unique], mean_expression[genes_sign_overlap]), 
  group = c(rep("all_genes", length(all_genes)), rep("sqtlseeker", length(genes_sign_sqtlseeker_unique)), rep("drimseq", length(genes_sign_drimseq_unique)), rep("overlap", length(genes_sign_overlap))))

ggdf$group <- factor(ggdf$group, levels = c("all_genes", "overlap", "sqtlseeker", "drimseq"))

ggdf <- ggdf[ggdf$mean_expression > 0, ,drop = FALSE]



ggp <- ggplot(ggdf, aes(x = log10(mean_expression), color = group, group = group)) +
  geom_density(size = 2) +
  theme_bw() +
  xlab("Log10 of mean gene expression ") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = as.character(c("grey", "black", colors[c("sqtlseeker", "drimseq")])))

pdf(paste0(out_dir, "sign_sqtls_mean_expr.pdf"))
print(ggp)
dev.off()


### calculate the number of expressed transcripts per gene 
nr_trans <- by(counts, factor(counts_raw$geneId), function(x){
  
  x <- as.matrix(x)
  
  sum(rowSums(x > 10, na.rm = TRUE) > 5, na.rm = TRUE)
  
}, simplify = FALSE)

nr_trans <- unlist(nr_trans)


ggdf <- data.frame(nr_trans = c(nr_trans[all_genes], nr_trans[genes_sign_sqtlseeker_unique], nr_trans[genes_sign_drimseq_unique], nr_trans[genes_sign_overlap]), 
  group = c(rep("all_genes", length(all_genes)), rep("sqtlseeker", length(genes_sign_sqtlseeker_unique)), rep("drimseq", length(genes_sign_drimseq_unique)), rep("overlap", length(genes_sign_overlap))))

ggdf$group <- factor(ggdf$group, levels = c("all_genes", "overlap", "sqtlseeker", "drimseq"))



ggp <- ggplot(ggdf, aes(x = nr_trans, color = group, group = group)) +
  geom_density(size = 2) +
  theme_bw() +
  xlab("Number of expressed transcripts") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = as.character(c("grey", "black", colors[c("sqtlseeker", "drimseq")])))

pdf(paste0(out_dir, "sign_sqtls_nr_trans.pdf"))
print(ggp)
dev.off()


############################################################################
# Check how many sQTLs is within exons
############################################################################

gtf <- import(path_gtf)

## keep exon regions for protein coding genes
keep <- mcols(gtf)$gene_type == "protein_coding" & mcols(gtf)$type == "exon"

gtf_exon <- gtf[keep, ]


freq_within_exons <- function(sqlt_list, gtf_exon, BPPARAM){
  
  gene_id <- strsplit2(sqlt_list, ":")[, 1]
  genes <- unique(gene_id)
  
  sqlt_per_gene <- split(sqlt_list, factor(gene_id, levels = genes))
  
  
  freq_list <- bplapply(1:length(sqlt_per_gene), function(i){
    # i = 1
    
    x <- sqlt_per_gene[[i]]
    
    name_split <- strsplit2(x, "_")
    start_snp <- as.numeric(name_split[, 3])
    
    snp_ranges <- GRanges(Rle(paste0("chr", name_split[, 2])), IRanges(start_snp, start_snp))
    gene_ranges <- gtf_exon[mcols(gtf_exon)$gene_id == genes[i], ]
    
    variantMatch <- GenomicRanges::findOverlaps(snp_ranges, gene_ranges, select = "first")
    
    return(!is.na(variantMatch))
    
  }, BPPARAM = BPPARAM)
  
  freq <- unlist(freq_list)
  freq <- mean(freq, na.rm = TRUE)
  
  return(freq)
  
}

non_sqtl <- setdiff(all_sqtls, union(sqtls_sign_sqtlseeker, sqtls_sign_drimseq))

freq_non_sqtl <- freq_within_exons(sqlt_list = non_sqtl, gtf_exon, BPPARAM = BPPARAM)
freq_sqtls_sign_overlap <- freq_within_exons(sqlt_list = sqtls_sign_overlap, gtf_exon, BPPARAM = BPPARAM)
freq_sqtls_sign_sqtlseeker_unique <- freq_within_exons(sqlt_list = sqtls_sign_sqtlseeker_unique, gtf_exon, BPPARAM = BPPARAM)
freq_sqtls_sign_drimseq_unique <- freq_within_exons(sqlt_list = sqtls_sign_drimseq_unique, gtf_exon, BPPARAM = BPPARAM)


freq_summary <- data.frame(set = c("non_sqtl", "overlap", "sqtlseeker", "drimseq"), freq_within_exon = c(freq_non_sqtl, freq_sqtls_sign_overlap, freq_sqtls_sign_sqtlseeker_unique, freq_sqtls_sign_drimseq_unique))

write.table(freq_summary, paste0(out_dir, "sign_sqtls_freq_within_exon.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


############################################################################
# Check how many sQTLs is within exons
############################################################################


dist_closest_exon <- function(sqlt_list, gtf_exon, BPPARAM){
  
  gene_id <- strsplit2(sqlt_list, ":")[, 1]
  genes <- unique(gene_id)
  
  sqlt_per_gene <- split(sqlt_list, factor(gene_id, levels = genes))
  
  
  dist_list <- bplapply(1:length(sqlt_per_gene), function(i){
    # i = 1
    
    x <- sqlt_per_gene[[i]]
    
    name_split <- strsplit2(x, "_")
    start_snp <- as.numeric(name_split[, 3])
    
    snp_ranges <- GRanges(Rle(paste0("chr", name_split[, 2])), IRanges(start_snp, start_snp))
    gene_ranges <- gtf_exon[mcols(gtf_exon)$gene_id == genes[i], ]
    
    variantMatch <- GenomicRanges::findOverlaps(snp_ranges, gene_ranges, select = "first")
    
    ### Set NA for snps that are within exons
    dist <- rep(NA, length(variantMatch))

    for(j in which(is.na(variantMatch))){
      # j = 1
      
      dist[j] <- min(abs(c(start_snp[j] - start(gene_ranges), start_snp[j] - end(gene_ranges))))
      
      }
    
    return(dist)

  }, BPPARAM = BPPARAM)
  
  dist <- unlist(dist_list)
  
  return(dist)
  
}


dist_non_sqtl <- dist_closest_exon(sqlt_list = non_sqtl, gtf_exon, BPPARAM = BPPARAM)
dist_sqtls_sign_overlap <- dist_closest_exon(sqlt_list = sqtls_sign_overlap, gtf_exon, BPPARAM = BPPARAM)
dist_sqtls_sign_sqtlseeker_unique <- dist_closest_exon(sqlt_list = sqtls_sign_sqtlseeker_unique, gtf_exon, BPPARAM = BPPARAM)
dist_sqtls_sign_drimseq_unique <- dist_closest_exon(sqlt_list = sqtls_sign_drimseq_unique, gtf_exon, BPPARAM = BPPARAM)


ggdf <- data.frame(dist = c(dist_non_sqtl, dist_sqtls_sign_overlap, dist_sqtls_sign_sqtlseeker_unique, dist_sqtls_sign_drimseq_unique), 
  group = c(rep("non_sqtl", length(dist_non_sqtl)), rep("overlap", length(dist_sqtls_sign_overlap)), rep("sqtlseeker", length(dist_sqtls_sign_sqtlseeker_unique)), rep("drimseq", length(dist_sqtls_sign_drimseq_unique))))

ggdf <- ggdf[complete.cases(ggdf), , drop = FALSE]

ggdf$group <- factor(ggdf$group, levels = c("non_sqtl", "overlap", "sqtlseeker", "drimseq"))



ggp <- ggplot(ggdf, aes(x = dist, color = group, group = group)) +
  stat_ecdf(size = 2) +
  theme_bw() +
  xlab("Distance to the closest exon") +
  ylab("Cumulative proportion") +
  coord_cartesian(xlim = c(0, 5000)) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = as.character(c("grey", "black", colors[c("sqtlseeker", "drimseq")])))

pdf(paste0(out_dir, "sign_sqtls_dist_closest_exon_cdf.pdf"))
print(ggp)
dev.off()










































