##############################################################################
# <<geuvadis_drimseq_0_3_3_comparison_permutations_adjsnps.R>>

# BioC 3.2
# Created 29 Feb 2016 
# Modified 21 Apr 2016

# For DRIMSeq merge results from all chromosomes; Calculate adjusted p-values for drimseq and sqtlseeker
# Plot histograms of p-values
# Plot venn diagrams
# Create a table with number of tested gene-snp pairs and significant sQTL per method

##############################################################################
Sys.time()
##############################################################################

library(DRIMSeq)
library(ggplot2)
library(iCOBRA)
library(plyr)
library(Gviz)
library(GenomicFeatures)
library(tools)
library(rtracklayer)
library(GenomicRanges)
library(limma)
library(matrixStats)

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/Shared/data/seq/geuvadis'
population='CEU'
method_out='drimseq_0_3_3_analysis_permutations_all_genes'
comparison_out='drimseq_0_3_3_comparison_permutations_all_genes_adjsnps'

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
print(population)
print(method_out)
print(comparison_out)

##############################################################################

setwd(rwd)

method_out <- paste0(method_out, "/")
comparison_out <- paste0(comparison_out, "/")

out_dir <- comparison_out
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


### colors

load(paste0(rwd, "/", comparison_out, "colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)

##############################################################################
# merge results 
##############################################################################

results <- list()

#####################################
### sqtlseeker results
#####################################


res_tmp <- read.table(paste0("sqtlseeker_2_1_analysis/results/", population, "_results_all.txt"), header = TRUE, as.is = TRUE)
head(res_tmp)

colnames(res_tmp) <- c("gene_id", "snp_id", "F", "nb.groups", "md", "tr.first", "tr.second", "nb.perms", "pvalue")
res_tmp <- unique(res_tmp)


res_tmp$gene_snp <- paste0(res_tmp$gene_id, ":", res_tmp$snp_id)

# res_tmp$adj_pvalue <- qvalue::qvalue(res_tmp$pvalue)$qvalues

res_tmp$adj_pvalue <- p.adjust(res_tmp$pvalue, method = "BH")


write.table(res_tmp, paste0(out_dir, "results_sqtlseeker.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


results[["sqtlseeker"]] <- res_tmp



# dim(results[["sqtlseeker"]])
# 
# table(duplicated(results[["sqtlseeker"]][, "gene_snp"]))



#####################################
### Read in DRIMSeq results + adjust p-values based on p-values from all the chromosomes
#####################################


res_list <- lapply(1:22, function(chr){
  # chr = 1
  
  res <- read.table(paste0(method_out, population, "_chr",chr, "_results.txt"), header = TRUE, sep = "\t", as.is = TRUE)
  
  return(res)
  
})


res <- rbind.fill(res_list)

res <- res[!is.na(res$pvalue), ]

res$gene_block <- paste0(res$gene_id, ":", res$block_id)
res$gene_snp <- paste0(res$gene_id, ":", res$snp_id)


### Redo the permutation adjustment of p-values based on p-values from all the genes

# Extract p-values from permutations
pval_perm_list <- lapply(1:22, function(chr){
  # chr = 1
  
  load(paste0(method_out, population, "_chr",chr, "_d.Rdata"))
  
  ### Permuted p-values are for blocks
  
  pval_perm_block <- data.frame(d@pvalues_permutations)
  pval_perm_block$block_id <- rownames(d@genotypes@unlistData)
  pval_perm_block$gene_id <- rep(names(d@genotypes), times = width(d@genotypes))
    
  ### Repeat results for blocks with multiple SNPs 

  pval_perm_block_spl <- split(pval_perm_block, factor(pval_perm_block$gene_id, levels = names(d@genotypes)))
  inds <- 1:length(pval_perm_block_spl)
  
  pval_perm_snp <- lapply(inds, function(i){
    # i = 1
    
    pval_perm_block_gene <- pval_perm_block_spl[[i]]
    
    blo <- d@blocks[[i]]
    matching <- match(blo[, "block_id"], pval_perm_block_gene[, "block_id"])
    
    pval_perm_snp_gene <- pval_perm_block_gene[matching, ]
    
    return(pval_perm_snp_gene)
    
  })
  
  pval_perm_snp <- do.call(rbind, pval_perm_snp)
  
  pval_perm_snp <- pval_perm_snp[, -grep("gene_id|block_id", colnames(pval_perm_snp))]
  
  dim(pval_perm_snp)
  dim(d@results)

  return(pval_perm_snp)
  
})



pval_perm <- do.call(rbind, pval_perm_list)

pval_perm <- as.matrix(pval_perm)

pval <- res$pvalue
nas <- is.na(pval)
pval <- pval[!nas]
pval <- factor(pval)


### Count how many pval_permuted is lower than pval from the model
pval_perm <- c(pval_perm)
nas_perm <- is.na(pval_perm)
pval_perm <- pval_perm[!nas_perm]
pval_perm_cut <- cut(pval_perm, c(-1, levels(pval), 2), right=FALSE)
pval_perm_sum <- table(pval_perm_cut)
pval_perm_cumsum <- cumsum(pval_perm_sum)[-length(pval_perm_sum)]
names(pval_perm_cumsum) <- levels(pval)
sum_sign_pval <- pval_perm_cumsum[pval]

nr_perm_tot <- length(pval_perm)
pval_adj <- (sum_sign_pval + 1) / (nr_perm_tot + 1)


res$pvalue_perm_new <- NA
res$pvalue_perm_new[!nas] <- pval_adj

res$adj_pvalue_perm_new <- p.adjust(res$pvalue_perm_new, method = "BH")


### Replace current p-values with one adjusted by permutations

res$pvalue <- res$pvalue_perm_new
res$adj_pvalue <- res$adj_pvalue_perm_new


### Remove the permutation columns

res <- res[, -grep("perm", colnames(res))]



write.table(res, paste0(out_dir, "results_drimseq.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


results[["drimseq"]] <- res





##########################################################################
### histograms of p-values
### scatterplot of p-values
##########################################################################


ggp <- DRIMSeq:::dm_plotPvalues(pvalues = results[["sqtlseeker"]][, "pvalue"]) +
  geom_histogram(breaks = seq(0,1, by = 0.01), fill = colors["sqtlseeker"])

pdf(paste0(out_dir, "sqtlseeker_hist_pvalues.pdf"))
print(ggp)
dev.off()


ggp <- DRIMSeq:::dm_plotPvalues(pvalues = results[["drimseq"]][, "pvalue"]) +
  geom_histogram(breaks = seq(0,1, by = 0.01), fill = colors["drimseq"])

pdf(paste0(out_dir, "drimseq_hist_pvalues.pdf"))
print(ggp)
dev.off()



### Plot scatter

results_pval <- list()

pval_tmp <- results[["sqtlseeker"]][, c("gene_snp", "pvalue")]
colnames(pval_tmp) <- c("gene_snp", "sqtlseeker")
results_pval[["sqtlseeker"]] <- pval_tmp


pval_tmp <- results[["drimseq"]][, c("gene_snp", "pvalue")]
colnames(pval_tmp) <- c("gene_snp", "drimseq")
results_pval[["drimseq"]] <- pval_tmp


results_pval <- Reduce(function(...) merge(..., by = "gene_snp", all=TRUE, sort = FALSE), results_pval)
rownames(results_pval) <- results_pval$gene_snp
results_pval <- results_pval[, -1]


results_pval_log <- -log10(results_pval)




pdf(paste0(out_dir, "scatter_pvalues.pdf"))
smoothScatter(results_pval[, "drimseq"], results_pval[, "sqtlseeker"], xlab = "p-value DRIMSeq", ylab = "p-value sQTLseekeR")
dev.off()

pdf(paste0(out_dir, "scatter_pvalues_log.pdf"))
smoothScatter(results_pval_log[, "drimseq"], results_pval_log[, "sqtlseeker"], xlab = "-log10(p-value) DRIMSeq", ylab = "-log10(p-value) sQTLseekeR")
dev.off()

png(paste0(out_dir, "scatter_pvalues_log.png"))
smoothScatter(results_pval_log[, "drimseq"], results_pval_log[, "sqtlseeker"], xlab = "-log10(p-value) DRIMSeq", ylab = "-log10(p-value) sQTLseekeR", nrpoints = nrow(results_pval_log))
dev.off()


# ggp <- ggplot(results_pval_log, aes(x = drimseq, y = sqtlseeker)) +
#   stat_binhex(bins = 100) 
# 
# pdf(paste0(out_dir, "binhex_pvalues.pdf"))
# print(ggp)
# dev.off()



# ggp <- ggplot(results_pval, aes(x = drimseq, y = sqtlseeker)) +
#   geom_point(alpha = 0.3) +
#   geom_abline(intercept = 0, slope = 1, color = "orange") +
#   xlab("p-value DRIMSeq") +
#   ylab("p-value sQTLseekeR")
# 
# png(paste0(out_dir, "scatter_pvalues.png"))
# print(ggp)
# dev.off()
# 
# 
# 
# ggp <- ggplot(results_pval_log, aes(x = drimseq, y = sqtlseeker)) +
#   geom_point(alpha = 0.3) +
#   geom_abline(intercept = 0, slope = 1, color = "orange") +
#   xlab("-log10(p-value) DRIMSeq") +
#   ylab("-log10(p-value) sQTLseekeR")
# 
# png(paste0(out_dir, "scatter_pvalues_log.png"))
# print(ggp)
# dev.off()


#######################################################
# merge results for iCOBRA - results_padj
#######################################################


results_padj <- list()

padj_tmp <- results[["sqtlseeker"]][, c("gene_snp", "adj_pvalue")]
colnames(padj_tmp) <- c("gene_snp", "sqtlseeker")
results_padj[["sqtlseeker"]] <- padj_tmp


padj_tmp <- results[["drimseq"]][, c("gene_snp", "adj_pvalue")]
colnames(padj_tmp) <- c("gene_snp", "drimseq")
results_padj[["drimseq"]] <- padj_tmp


results_padj <- Reduce(function(...) merge(..., by = "gene_snp", all=TRUE, sort = FALSE), results_padj)
rownames(results_padj) <- results_padj$gene_snp


results_padj <- results_padj[, colnames(results_padj) %in% colors_df$methods, drop = FALSE]

keep_methods <- colors_df$methods %in% colnames(results_padj)

colors <- colors[keep_methods]
colors_df <- colors_df[keep_methods, , drop = FALSE]




##########################################################################
###  use iCOBRA for overlaps
##########################################################################


summary <- data.frame(method = c("drimseq", "sqtlseeker", "overlap"))




# Significant gene-SNPs


cobradata <- COBRAData(padj = results_padj)

cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 0.05)

basemethods(cobraperf)

colorscheme <- colors[basemethods(cobraperf)]

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colorscheme, incltruth = FALSE)


pdf(paste0(out_dir, "/venn_gene_snp_sign.pdf"))
plot_overlap(cobraplot, cex=c(1.2,1,0.7))
dev.off()


pdf(paste0(out_dir, "/upset_gene_snp_sign.pdf"))
plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 15)
dev.off()



overlap <- cobraplot@overlap

summary$gene_snp_sign <- c(colSums(overlap), sum(rowSums(overlap == 1) == 2))



### Significant genes

# Keep minimal p-value per gene
gene_ids <- factor(strsplit2(rownames(results_padj), ":")[, 1])

results_padj_split <- by(results_padj, gene_ids, function(x){
  
  colMins(as.matrix(x), na.rm = TRUE)
  
}, simplify = FALSE)


results_padj_gene <- data.frame(do.call(rbind, results_padj_split))
colnames(results_padj_gene) <- colnames(results_padj)

# Use iCOBRA
cobradata <- COBRAData(padj = results_padj_gene)

cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 0.05)

basemethods(cobraperf)

colorscheme <- colors[basemethods(cobraperf)]

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colorscheme, incltruth = FALSE)


pdf(paste0(out_dir, "/venn_gene_sign.pdf"))
plot_overlap(cobraplot, cex=c(1.2,1,0.7))
dev.off()


pdf(paste0(out_dir, "/upset_gene_sign.pdf"))
plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 15)
dev.off()



overlap <- cobraplot@overlap

summary$gene_sign <- c(colSums(overlap), sum(rowSums(overlap == 1) == 2))



# All tested gene-SNPs

cobradata <- COBRAData(padj = results_padj)

cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 1.1)

basemethods(cobraperf)


colorscheme <- colors[basemethods(cobraperf)]

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colorscheme, incltruth = FALSE)


pdf(paste0(out_dir, "/venn_gene_snp_all.pdf"))
plot_overlap(cobraplot, cex=c(1.2,1,0.7))
dev.off()


pdf(paste0(out_dir, "/upset_gene_snp_all.pdf"))
plot_upset(cobraplot, order.by = "degree", empty.intersections = "on")
dev.off()


overlap <- cobraplot@overlap

summary$gene_snp_all <- c(colSums(overlap), sum(rowSums(overlap == 1) == 2))


# All tested genes


cobradata <- COBRAData(padj = results_padj_gene)

cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 1.1)

basemethods(cobraperf)


colorscheme <- colors[basemethods(cobraperf)]

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colorscheme, incltruth = FALSE)


pdf(paste0(out_dir, "/venn_gene_all.pdf"))
plot_overlap(cobraplot, cex=c(1.2,1,0.7))
dev.off()


pdf(paste0(out_dir, "/upset_gene_all.pdf"))
plot_upset(cobraplot, order.by = "degree", empty.intersections = "on")
dev.off()


overlap <- cobraplot@overlap

summary$gene_all <- c(colSums(overlap), sum(rowSums(overlap == 1) == 2))




write.table(summary, paste0(out_dir, "/summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
















































