##############################################################################
# <<geuvadis_drimseq_0_3_3_comparison_run.R>>

# BioC 3.2
# Created 29 Feb 2016 
# Modified 9 Apr 2016

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

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/Shared/data/seq/geuvadis'
population='CEU'
method_out='drimseq_0_3_3_analysis_permutations'
comparison_out='drimseq_0_3_3_comparison_permutations'

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


res_uniq <- res[!duplicated(res$gene_block), ]

res_uniq$adj_pvalue <- p.adjust(res_uniq$pvalue, method = "BH")

mm <- match(res$gene_block, res_uniq$gene_block)
res$adj_pvalue <- res_uniq$adj_pvalue[mm]


results[["drimseq"]] <- res





##########################################################################
### histograms of p-values
### scatterplot of p-values
##########################################################################


ggp <- DRIMSeq:::dm_plotPvalues(pvalues = results[["sqtlseeker"]][, "pvalue"]) +
  geom_histogram(breaks = seq(0,1, by = 0.01), fill = "deeppink4")

pdf(paste0(out_dir, "sqtlseeker_hist_pvalues.pdf"))
print(ggp)
dev.off()


ggp <- DRIMSeq:::dm_plotPvalues(pvalues = results[["drimseq"]][, "pvalue"]) +
  geom_histogram(breaks = seq(0,1, by = 0.01), fill = "deeppink4")

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


# All tested cases

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

overlap <- cobraplot@overlap[!duplicated(strsplit2(rownames(cobraplot@overlap), ":")[, 1]), ]

summary$gene_all <- c(colSums(overlap), sum(rowSums(overlap == 1) == 2))


pdf(paste0(out_dir, "/venn_gene_all.pdf"))
vennDiagram(overlap, circle.col = colorscheme)
dev.off()



# Significant cases


cobradata <- COBRAData(padj = results_padj)

cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 0.05)

basemethods(cobraperf)

colorscheme <- colors[basemethods(cobraperf)]

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colorscheme, incltruth = FALSE)


pdf(paste0(out_dir, "/venn_gene_snp_sign.pdf"))
plot_overlap(cobraplot, cex=c(1.2,1,0.7))
dev.off()


pdf(paste0(out_dir, "/upset_gene_snp_sign.pdf"))
plot_upset(cobraplot, order.by = "degree", empty.intersections = "on")
dev.off()



overlap <- cobraplot@overlap

summary$gene_snp_sign <- c(colSums(overlap), sum(rowSums(overlap == 1) == 2))


overlap_split <- split.data.frame(overlap, strsplit2(rownames(cobraplot@overlap), ":")[, 1])

overlap_list <- lapply(overlap_split, function(i){
  
  colSums(i) > 0
  
})



overlap <- do.call(rbind, overlap_list)

summary$gene_sign <- c(colSums(overlap), sum(rowSums(overlap == 1) == 2))

pdf(paste0(out_dir, "/venn_gene_sign.pdf"))
vennDiagram(overlap, circle.col = colorscheme)
dev.off()




write.table(summary, paste0(out_dir, "/summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
















































