##############################################################################
## <<geuvadis_drimseq_0_3_3_biol_comparison_plots.R>>

# BioC 3.2
# Created 21 Apr 2016
# Modified 17 Oct 2016


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
# comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_drimseq_counts'
# FDR=0.05
# path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf'
# workers=5

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

strip_text_size <- 16
text_size <- 18

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


res <- read.table(paste0(comparison_out, population, "_results_sqtlseeker.txt"), header = TRUE, as.is = TRUE)
head(res)

res <- res[!is.na(res$adj_pvalue), , drop = FALSE]

results[["sqtlseeker"]] <- res
metadata[["sqtlseeker"]] <- data.frame(method_name = "sqtlseeker", stringsAsFactors = FALSE)


#####################################
### DRIMSeq results 
#####################################

res <- read.table(paste0(comparison_out, population, "_results_drimseq.txt"), header = TRUE, as.is = TRUE)
head(res)

res <- res[!is.na(res$adj_pvalue), , drop = FALSE]

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



### Use unique genes per method
ggdf <- data.frame(mean_expression = c(mean_expression[all_genes], mean_expression[genes_sign_sqtlseeker_unique], mean_expression[genes_sign_drimseq_unique], mean_expression[genes_sign_overlap]),
  group = c(rep("all_genes", length(all_genes)), rep("sqtlseeker_unique", length(genes_sign_sqtlseeker_unique)), rep("drimseq_unique", length(genes_sign_drimseq_unique)), rep("overlap", length(genes_sign_overlap))))

ggdf$group <- factor(ggdf$group, levels = c("all_genes", "overlap", "sqtlseeker_unique", "drimseq_unique"), labels = paste0(c("all_genes", "overlap", "sqtlseeker_unique", "drimseq_unique"), " (", c(length(all_genes), length(genes_sign_overlap), length(genes_sign_sqtlseeker_unique), length(genes_sign_drimseq_unique)), ")"))

ggdf <- ggdf[ggdf$mean_expression > 0, ,drop = FALSE]


ggp <- ggplot(ggdf, aes(x = log10(mean_expression), color = group, group = group)) +
  geom_density(size = 2) +
  theme_bw() +
  xlab("Log10 of mean gene expression ") +
  theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = text_size), legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = as.character(c("grey", "gray50", colors[c("sqtlseeker", "drimseq")]))) +
  guides(color = guide_legend(nrow = 2))

pdf(paste0(out_dir, population, "_sign_sqtls_mean_expr.pdf"))
print(ggp)
dev.off()




### Use all genes per method
ggdf <- data.frame(mean_expression = c(mean_expression[all_genes], mean_expression[genes_sign_sqtlseeker], mean_expression[genes_sign_drimseq], mean_expression[genes_sign_overlap]),
  group = c(rep("all_genes", length(all_genes)), rep("sqtlseeker", length(genes_sign_sqtlseeker)), rep("drimseq", length(genes_sign_drimseq)), rep("overlap", length(genes_sign_overlap))))

ggdf$group <- factor(ggdf$group, levels = c("all_genes", "overlap", "sqtlseeker", "drimseq"), labels = paste0(c("all_genes", "overlap", "sqtlseeker", "drimseq"), " (", c(length(all_genes), length(genes_sign_overlap), length(genes_sign_sqtlseeker), length(genes_sign_drimseq)), ")"))

ggdf <- ggdf[ggdf$mean_expression > 0, ,drop = FALSE]



ggp <- ggplot(ggdf, aes(x = log10(mean_expression), color = group, group = group)) +
  geom_density(size = 2) +
  theme_bw() +
  xlab("Log10 of mean gene expression ") +
  theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = text_size), legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = as.character(c("grey", "gray50", colors[c("sqtlseeker", "drimseq")]))) +
  guides(color = guide_legend(nrow = 2))

pdf(paste0(out_dir, population, "_sign_sqtls_mean_expr2.pdf"))
print(ggp)
dev.off()




### calculate the number of expressed transcripts per gene
nr_trans <- by(counts, factor(counts_raw$geneId), function(x){

  x <- as.matrix(x)

  sum(rowSums(x > 10, na.rm = TRUE) > 5, na.rm = TRUE)

}, simplify = FALSE)

nr_trans <- unlist(nr_trans)



### Use unique genes per method
ggdf <- data.frame(nr_trans = c(nr_trans[all_genes], nr_trans[genes_sign_sqtlseeker_unique], nr_trans[genes_sign_drimseq_unique], nr_trans[genes_sign_overlap]),
  group = c(rep("all_genes", length(all_genes)), rep("sqtlseeker_unique", length(genes_sign_sqtlseeker_unique)), rep("drimseq_unique", length(genes_sign_drimseq_unique)), rep("overlap", length(genes_sign_overlap))))

ggdf$group <- factor(ggdf$group, levels = c("all_genes", "overlap", "sqtlseeker_unique", "drimseq_unique"), labels = paste0(c("all_genes", "overlap", "sqtlseeker_unique", "drimseq_unique"), " (", c(length(all_genes), length(genes_sign_overlap), length(genes_sign_sqtlseeker_unique), length(genes_sign_drimseq_unique)), ")"))


ggp <- ggplot(ggdf, aes(x = nr_trans, color = group, group = group)) +
  geom_density(size = 2) +
  theme_bw() +
  xlab("Number of expressed transcripts") +
  theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = text_size), legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = as.character(c("grey", "gray50", colors[c("sqtlseeker", "drimseq")]))) +
  guides(color = guide_legend(nrow = 2))

pdf(paste0(out_dir, population, "_sign_sqtls_nr_trans.pdf"))
print(ggp)
dev.off()



### Use all genes per method
ggdf <- data.frame(nr_trans = c(nr_trans[all_genes], nr_trans[genes_sign_sqtlseeker], nr_trans[genes_sign_drimseq], nr_trans[genes_sign_overlap]),
  group = c(rep("all_genes", length(all_genes)), rep("sqtlseeker", length(genes_sign_sqtlseeker)), rep("drimseq", length(genes_sign_drimseq)), rep("overlap", length(genes_sign_overlap))))

ggdf$group <- factor(ggdf$group, levels = c("all_genes", "overlap", "sqtlseeker", "drimseq"), labels = paste0(c("all_genes", "overlap", "sqtlseeker", "drimseq"), " (", c(length(all_genes), length(genes_sign_overlap), length(genes_sign_sqtlseeker), length(genes_sign_drimseq)), ")"))


ggp <- ggplot(ggdf, aes(x = nr_trans, color = group, group = group)) +
  geom_density(size = 2) +
  theme_bw() +
  xlab("Number of expressed transcripts") +
  theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = text_size), legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = as.character(c("grey", "gray50", colors[c("sqtlseeker", "drimseq")]))) +
  guides(color = guide_legend(nrow = 2))

pdf(paste0(out_dir, population, "_sign_sqtls_nr_trans2.pdf"))
print(ggp)
dev.off()




### Plot a scatter of nr of transcripts versus mean gene expression

ggdf <- data.frame(mean_expression = c(mean_expression[genes_sign_sqtlseeker_unique], mean_expression[genes_sign_drimseq_unique], mean_expression[genes_sign_overlap]),
  nr_trans = c(nr_trans[genes_sign_sqtlseeker_unique], nr_trans[genes_sign_drimseq_unique], nr_trans[genes_sign_overlap]),
  group = c(rep("sqtlseeker_unique", length(genes_sign_sqtlseeker_unique)), rep("drimseq_unique", length(genes_sign_drimseq_unique)), rep("overlap", length(genes_sign_overlap))))

ggdf$group <- factor(ggdf$group, levels = c("overlap", "sqtlseeker_unique", "drimseq_unique"), labels = paste0(c("overlap", "sqtlseeker_unique", "drimseq_unique"), " (", c(length(genes_sign_overlap), length(genes_sign_sqtlseeker_unique), length(genes_sign_drimseq_unique)), ")"))


ggp <- ggplot(ggdf, aes(x = log10(mean_expression), y = nr_trans, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  xlab("Log10 of mean gene expression ") +
  ylab("Number of expressed transcripts") +
  theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = text_size), legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = as.character(c("gray50", colors[c("sqtlseeker", "drimseq")]))) +
  guides(color = guide_legend(nrow = 2))

pdf(paste0(out_dir, population, "_sign_sqtls_scatter.pdf"))
print(ggp)
dev.off()



############################################################################
# Check how many sQTLs are within exons
############################################################################

gtf <- import(path_gtf)

## keep exon regions for protein coding genes
keep <- mcols(gtf)$gene_type == "protein_coding" & mcols(gtf)$type == "exon"

gtf_exon <- gtf[keep, ]



freq_within_ranges <- function(sqlt_list, gtf, BPPARAM){
  
  gene_id <- strsplit2(sqlt_list, ":")[, 1]
  genes <- unique(gene_id)
  
  sqlt_per_gene <- split(sqlt_list, factor(gene_id, levels = genes))
  
  
  freq_list <- bplapply(1:length(sqlt_per_gene), function(i){
    # i = 1
    
    x <- sqlt_per_gene[[i]]
    
    name_split <- strsplit2(x, "_")
    start_snp <- as.numeric(name_split[, 3])
    
    snp_ranges <- GRanges(Rle(paste0("chr", name_split[, 2])), IRanges(start_snp, start_snp))
    ranges <- gtf[mcols(gtf)$gene_id == genes[i], ]
    
    variantMatch <- GenomicRanges::findOverlaps(snp_ranges, ranges, select = "first")
    
    return(!is.na(variantMatch))
    
  }, BPPARAM = BPPARAM)
  
  freq <- unlist(freq_list)
  freq <- mean(freq, na.rm = TRUE)
  
  return(freq)
  
}

non_sqtl <- setdiff(all_sqtls, union(sqtls_sign_sqtlseeker, sqtls_sign_drimseq))

freq_non_sqtl <- freq_within_ranges(sqlt_list = non_sqtl, gtf_exon, BPPARAM = BPPARAM)
freq_sqtls_sign_overlap <- freq_within_ranges(sqlt_list = sqtls_sign_overlap, gtf_exon, BPPARAM = BPPARAM)

freq_sqtls_sign_sqtlseeker_unique <- freq_within_ranges(sqlt_list = sqtls_sign_sqtlseeker_unique, gtf_exon, BPPARAM = BPPARAM)
freq_sqtls_sign_drimseq_unique <- freq_within_ranges(sqlt_list = sqtls_sign_drimseq_unique, gtf_exon, BPPARAM = BPPARAM)

freq_sqtls_sign_sqtlseeker <- freq_within_ranges(sqlt_list = sqtls_sign_sqtlseeker, gtf_exon, BPPARAM = BPPARAM)
freq_sqtls_sign_drimseq <- freq_within_ranges(sqlt_list = sqtls_sign_drimseq, gtf_exon, BPPARAM = BPPARAM)


freq_summary <- data.frame(set = c("non_sqtl", "overlap", "sqtlseeker", "drimseq", "sqtlseeker_unique", "drimseq_unique"), freq_within_exon = c(freq_non_sqtl, freq_sqtls_sign_overlap, freq_sqtls_sign_sqtlseeker, freq_sqtls_sign_drimseq, freq_sqtls_sign_sqtlseeker_unique, freq_sqtls_sign_drimseq_unique))


write.table(freq_summary, paste0(out_dir, population, "_sign_sqtls_freq_within_exon.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


############################################################################
# Check how many sQTLs are within splice sites 
############################################################################


mcols(gtf_exon)$transcript_id <- factor(mcols(gtf_exon)$transcript_id)
transcripts <- levels(mcols(gtf_exon)$transcript_id)


## get the splice sites per transcript
gtf_sss <- bplapply(1:nlevels(mcols(gtf_exon)$transcript_id), function(i){
  # i = 1
  gtf_transcript <- gtf_exon[mcols(gtf_exon)$transcript_id == transcripts[i], ]
  
  # have to give start, otherwise the first range is from 1
  introns <- gaps(gtf_transcript, start = min(start(gtf_transcript)))
  
  ## Extract the splice stites
  donor <- introns
  end(donor) <- start(donor) + 1
  acceptor <- introns
  start(acceptor) <- end(acceptor) - 1
  
  splice_sites <- c(donor, acceptor)
  mcols(splice_sites) <- data.frame(transcript_id = transcripts[i], gene_id = mcols(gtf_transcript)$gene_id[1], stringsAsFactors = FALSE)
  
  return(splice_sites)
  
}, BPPARAM = BPPARAM)

gtf_sss <- do.call(c, gtf_sss)


non_sqtl <- setdiff(all_sqtls, union(sqtls_sign_sqtlseeker, sqtls_sign_drimseq))

freq_non_sqtl <- freq_within_ranges(sqlt_list = non_sqtl, gtf_sss, BPPARAM = BPPARAM)
freq_sqtls_sign_overlap <- freq_within_ranges(sqlt_list = sqtls_sign_overlap, gtf_sss, BPPARAM = BPPARAM)

freq_sqtls_sign_sqtlseeker_unique <- freq_within_ranges(sqlt_list = sqtls_sign_sqtlseeker_unique, gtf_sss, BPPARAM = BPPARAM)
freq_sqtls_sign_drimseq_unique <- freq_within_ranges(sqlt_list = sqtls_sign_drimseq_unique, gtf_sss, BPPARAM = BPPARAM)

freq_sqtls_sign_sqtlseeker <- freq_within_ranges(sqlt_list = sqtls_sign_sqtlseeker, gtf_sss, BPPARAM = BPPARAM)
freq_sqtls_sign_drimseq <- freq_within_ranges(sqlt_list = sqtls_sign_drimseq, gtf_sss, BPPARAM = BPPARAM)


freq_summary <- data.frame(set = c("non_sqtl", "overlap", "sqtlseeker", "drimseq", "sqtlseeker_unique", "drimseq_unique"), freq_within_exon = c(freq_non_sqtl, freq_sqtls_sign_overlap, freq_sqtls_sign_sqtlseeker, freq_sqtls_sign_drimseq, freq_sqtls_sign_sqtlseeker_unique, freq_sqtls_sign_drimseq_unique))


write.table(freq_summary, paste0(out_dir, population, "_sign_sqtls_freq_within_splice_sites.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



############################################################################
# Check the distance of nonexonic sQTLs to the closest exon
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

dist_sqtls_sign_sqtlseeker <- dist_closest_exon(sqlt_list = sqtls_sign_sqtlseeker, gtf_exon, BPPARAM = BPPARAM)
dist_sqtls_sign_drimseq <- dist_closest_exon(sqlt_list = sqtls_sign_drimseq, gtf_exon, BPPARAM = BPPARAM)




### Use unique genes per method
ggdf <- data.frame(dist = c(dist_non_sqtl, dist_sqtls_sign_overlap, dist_sqtls_sign_sqtlseeker_unique, dist_sqtls_sign_drimseq_unique),
  group = c(rep("non_sqtl", length(dist_non_sqtl)), rep("overlap", length(dist_sqtls_sign_overlap)), rep("sqtlseeker_unique", length(dist_sqtls_sign_sqtlseeker_unique)), rep("drimseq_unique", length(dist_sqtls_sign_drimseq_unique))))

ggdf <- ggdf[complete.cases(ggdf), , drop = FALSE]

ggdf$group <- factor(ggdf$group, levels = c("non_sqtl", "overlap", "sqtlseeker_unique", "drimseq_unique"), labels = paste0(c("non_sqtl", "overlap", "sqtlseeker_unique", "drimseq_unique"), " (", c(length(dist_non_sqtl), length(dist_sqtls_sign_overlap), length(dist_sqtls_sign_sqtlseeker_unique), length(dist_sqtls_sign_drimseq_unique)), ")"))


ggp <- ggplot(ggdf, aes(x = dist, color = group, group = group)) +
  stat_ecdf(size = 2) +
  theme_bw() +
  xlab("Distance to the closest exon") +
  ylab("Cumulative proportion") +
  coord_cartesian(xlim = c(0, 5000)) +
  theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = text_size), legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = as.character(c("grey", "gray50", colors[c("sqtlseeker", "drimseq")]))) +
  guides(color = guide_legend(nrow = 2))

pdf(paste0(out_dir, population, "_sign_sqtls_dist_closest_exon_cdf.pdf"))
print(ggp)
dev.off()



### Use all genes per method
ggdf <- data.frame(dist = c(dist_non_sqtl, dist_sqtls_sign_overlap, dist_sqtls_sign_sqtlseeker, dist_sqtls_sign_drimseq),
  group = c(rep("non_sqtl", length(dist_non_sqtl)), rep("overlap", length(dist_sqtls_sign_overlap)), rep("sqtlseeker", length(dist_sqtls_sign_sqtlseeker)), rep("drimseq", length(dist_sqtls_sign_drimseq))))

ggdf <- ggdf[complete.cases(ggdf), , drop = FALSE]

ggdf$group <- factor(ggdf$group, levels = c("non_sqtl", "overlap", "sqtlseeker", "drimseq"), labels = paste0(c("non_sqtl", "overlap", "sqtlseeker", "drimseq"), " (", c(length(dist_non_sqtl), length(dist_sqtls_sign_overlap), length(dist_sqtls_sign_sqtlseeker), length(dist_sqtls_sign_drimseq)), ")"))


ggp <- ggplot(ggdf, aes(x = dist, color = group, group = group)) +
  stat_ecdf(size = 2) +
  theme_bw() +
  xlab("Distance to the closest exon") +
  ylab("Cumulative proportion") +
  coord_cartesian(xlim = c(0, 5000)) +
  theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = text_size), legend.title = element_blank(), legend.position = "bottom") +
  scale_color_manual(values = as.character(c("grey", "gray50", colors[c("sqtlseeker", "drimseq")]))) +
  guides(color = guide_legend(nrow = 2))

pdf(paste0(out_dir, population, "_sign_sqtls_dist_closest_exon_cdf2.pdf"))
print(ggp)
dev.off()



############################################################################
# Check how many sQTLs is within 1kb GWAS
############################################################################

### Read in the GWAS data
path_gwas <- "data/gwas/gwas_catalog_v1.0-associations_e84_r2016-05-01.tsv"

gwas <- read.table(path_gwas, header = TRUE, sep = "\t", as.is = TRUE, quote = "", fill = TRUE)

### Read in the snp ID converter
path_id_convert <- "data/metadata/snp_id_convert.Rdata"
load(path_id_convert)


match_gwas <- match(gwas$SNPS, snp_id_convert[, 1])

gwas$snp_id <- snp_id_convert[match_gwas, 2]


### Keep only SNP (no indels, no NAs)

gwas_sub <- gwas[grep("snp_", gwas$snp_id), c("SNPS", "snp_id")]


### Create ranges around GWAS

gwas_split <- data.frame(strsplit2(gwas_sub$snp_id, "_"), stringsAsFactors = FALSE)

gwas_split[, 3] <- as.numeric(gwas_split[, 3])

gwas_ranges <- GRanges(Rle(paste0("chr", gwas_split[, 2])), IRanges(gwas_split[, 3], gwas_split[, 3]))

window <- 1001

gwas_ranges <- resize(gwas_ranges, window, fix = "center")



freq_within_gwas <- function(sqlt_list, gwas_ranges){

  snp_id <- unique(strsplit2(sqlt_list, ":")[, 2])

  snp_split <- data.frame(strsplit2(snp_id, "_"), stringsAsFactors = FALSE)

  start_snp <- as.numeric(snp_split[, 3])

  snp_ranges <- GRanges(Rle(paste0("chr", snp_split[, 2])), IRanges(start_snp, start_snp))

  variantMatch <- findOverlaps(snp_ranges, gwas_ranges, select = "first")

  freq <- !is.na(variantMatch)

  freq <- mean(freq, na.rm = TRUE)

  return(freq)

}

non_sqtl <- setdiff(all_sqtls, union(sqtls_sign_sqtlseeker, sqtls_sign_drimseq))

freq_non_sqtl <- freq_within_gwas(sqlt_list = non_sqtl, gwas_ranges)
freq_sqtls_sign_overlap <- freq_within_gwas(sqlt_list = sqtls_sign_overlap, gwas_ranges)

freq_sqtls_sign_sqtlseeker_unique <- freq_within_gwas(sqlt_list = sqtls_sign_sqtlseeker_unique, gwas_ranges)
freq_sqtls_sign_drimseq_unique <- freq_within_gwas(sqlt_list = sqtls_sign_drimseq_unique, gwas_ranges)

freq_sqtls_sign_sqtlseeker <- freq_within_gwas(sqlt_list = sqtls_sign_sqtlseeker, gwas_ranges)
freq_sqtls_sign_drimseq <- freq_within_gwas(sqlt_list = sqtls_sign_drimseq, gwas_ranges)



freq_summary <- data.frame(set = c("non_sqtl", "overlap", "sqtlseeker", "drimseq", "sqtlseeker_unique", "drimseq_unique"), freq_within_gwas = c(freq_non_sqtl, freq_sqtls_sign_overlap, freq_sqtls_sign_sqtlseeker, freq_sqtls_sign_drimseq, freq_sqtls_sign_sqtlseeker_unique, freq_sqtls_sign_drimseq_unique))

write.table(freq_summary, paste0(out_dir, population, "_sign_sqtls_freq_within_gwas.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)






































