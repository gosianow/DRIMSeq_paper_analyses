##############################################################################
# <<geuvadis_drimseq_0_3_3_comparison_permutations.R>>

# BioC 3.2
# Created 29 Feb 2016 
# Modified 12 Apr 2016

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
method_out='drimseq_0_3_3_analysis_permutations_all_genes'
comparison_out='drimseq_0_3_3_comparison_permutations_all_genes'
sqtlseeker_results='sqtlseeker_2_1_analysis'
FDR=0.05

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


res <- read.table(paste0(sqtlseeker_results, "/results/", population, "_results_all.txt"), header = TRUE, as.is = TRUE)
head(res)

colnames(res) <- c("gene_id", "snp_id", "F", "nb.groups", "md", "tr.first", "tr.second", "nb.perms", "pvalue")
res <- unique(res)


res$gene_snp <- paste0(res$gene_id, ":", res$snp_id)

# res$adj_pvalue <- qvalue::qvalue(res$pvalue)$qvalues

res$adj_pvalue <- p.adjust(res$pvalue, method = "BH")


write.table(res, paste0(out_dir, "results_sqtlseeker.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


results[["sqtlseeker"]] <- res



# dim(results[["sqtlseeker"]])
# 
# table(duplicated(results[["sqtlseeker"]][, "gene_snp"]))


### See how many genes were tested

length(unique(res$gene_id))


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

### Keep data for unique blocks
res_uniq <- res[!duplicated(res$gene_block), ]

### Redo the permutation adjustment of p-values based on p-values from all the genes


### Extract p-values from permutations
pval_perm_list <- lapply(1:22, function(chr){
  # chr = 1
  
  load(paste0(method_out, population, "_chr",chr, "_d.Rdata"))
  
  return(d@pvalues_permutations)
  
})

pval_perm <- do.call(rbind, pval_perm_list)


pval <- res_uniq$pvalue
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


res_uniq$pvalue_perm_new <- NA
res_uniq$pvalue_perm_new[!nas] <- pval_adj

res_uniq$adj_pvalue_perm_new <- p.adjust(res_uniq$pvalue_perm_new, method = "BH")


# pdf(paste0(out_dir, "pvalues_perm_new.pdf"))
# smoothScatter(res_uniq$pvalue_perm, res_uniq$pvalue_perm_new)
# dev.off()



### Remove the permutation columns

res <- res[, -grep("perm", colnames(res))]


### Replace current p-values with one adjusted by permutations

mm <- match(res$gene_block, res_uniq$gene_block)
res$pvalue <- res_uniq$pvalue_perm_new[mm]
res$adj_pvalue <- res_uniq$adj_pvalue_perm_new[mm]


write.table(res, paste0(out_dir, "results_drimseq.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


results[["drimseq"]] <- res


##########################################################################
### histograms of number of transcripts per gene (for genes that were returned in the results)
##########################################################################

### sqtlseeker

data <- read.table(paste0(sqtlseeker_results, "/data/trExpCount_", population, "_sqtlseeker_ratios.tsv"), header = TRUE, as.is = TRUE)


data <- data[data$geneId %in% results[["sqtlseeker"]]$gene_id, , drop = FALSE]


nr_transcripts_sqtlseeker <- as.numeric(table(data$geneId))

df <- data.frame(tt = nr_transcripts_sqtlseeker)

ggp <- ggplot(df, aes_string(x = "tt")) +
  geom_histogram(fill = colors["sqtlseeker"], breaks = seq(0, max(df$tt), by = 1)) +
  geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(length(df$tt), " genes   \n ", sum(df$tt) , " features   ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, size = 6) +
  xlab("Number of features per gene") +
  ylab("Frequency") +
  theme_bw() +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold")) +
  coord_cartesian(xlim = c(0, max(df$tt) + 2)) 
 

pdf(paste0(out_dir, "sqtlseeker_hist_features.pdf"))
print(ggp)
dev.off()



### drimseq


nr_transcripts <- lapply(1:22, function(chr){
  # chr = 1
  
  load(paste0(method_out, population, "_chr",chr, "_d.Rdata"))
  
  tt <- width(d@counts)
  
  return(tt)
  
})

nr_transcripts_drimseq <-  unlist(nr_transcripts)



df <- data.frame(tt = nr_transcripts_drimseq)

ggp <- ggplot(df, aes_string(x = "tt")) +
  geom_histogram(fill = colors["drimseq"], breaks = seq(0, max(df$tt), by = 1)) +
  geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(length(df$tt), " genes   \n ", sum(df$tt) , " features   ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, size = 6) +
  xlab("Number of features per gene") +
  ylab("Frequency") +
  theme_bw() +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold")) +
  coord_cartesian(xlim = c(0, max(df$tt) + 2)) 


pdf(paste0(out_dir, "drimseq_hist_features.pdf"))
print(ggp)
dev.off()



### both at one plot


df <- data.frame(tt = c(nr_transcripts_drimseq, nr_transcripts_sqtlseeker), method = c(rep("drimseq", length(nr_transcripts_drimseq)), rep("sqtlseeker", length(nr_transcripts_sqtlseeker))))


ggp <- ggplot(df, aes(x = tt, fill = method)) +
  geom_histogram(breaks = seq(0, max(df$tt), by = 1), alpha = 0.5, position="identity") +
  xlab("Number of features per gene") +
  ylab("Frequency") +
  theme_bw() +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold"),  legend.title = element_blank(), legend.position = "bottom") +
  scale_fill_manual(values = colors[order(colors, decreasing = TRUE)])

pdf(paste0(out_dir, "hist_features.pdf"))
print(ggp)
dev.off()


##########################################################################
### dispersion versus mean plot for drimseq
##########################################################################

dips_mean <- lapply(1:22, function(chr){
  # chr = 1
  
  load(paste0(method_out, population, "_chr",chr, "_d.Rdata"))
  
  w <- sapply(d@genewise_dispersion, length)
  
  mean_expression <- rep(d@mean_expression, w)
  nr_features <- rep(width(d@counts), w)
  
  genewise_dispersion <- unlist(d@genewise_dispersion)
  
  return(data.frame(mean_expression = mean_expression, genewise_dispersion = genewise_dispersion, nr_features = nr_features))
  
})

dips_mean <- rbind.fill(dips_mean)

dips_mean <- unique(dips_mean)

ggp <- DRIMSeq:::dm_plotDispersion(genewise_dispersion = dips_mean$genewise_dispersion, mean_expression = dips_mean$mean_expression, nr_features = dips_mean$nr_features, common_dispersion = NULL)


pdf(paste0(out_dir, "drimseq_disversion_versus_mean.pdf"))
print(ggp)
dev.off()



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


### both

df <- data.frame(pvalues = c(results[["drimseq"]][, "pvalue"], results[["sqtlseeker"]][, "pvalue"]), method = c(rep("drimseq", length(results[["drimseq"]][, "pvalue"])), rep("sqtlseeker", length(results[["sqtlseeker"]][, "pvalue"]))))


ggp <- ggplot(df, aes(x = pvalues, fill = method)) +
  theme_bw() +
  xlab("p-values") +
  ylab("Frequency") +
  geom_histogram(breaks = seq(0, 1, by = 0.01), alpha = 0.5, position="identity") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=16, face="bold"), legend.title = element_blank(), legend.position = "bottom") +
  coord_cartesian(xlim = c(0, 1)) +
  scale_fill_manual(values = colors[order(colors, decreasing = TRUE)])

pdf(paste0(out_dir, "hist_pvalues.pdf"))
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


results_padj <- results_padj[, colors_df$methods, drop = FALSE]


##########################################################################
###  use iCOBRA for overlaps
##########################################################################


summary <- data.frame(method = c(colors_df$methods, "overlap"))



# Significant gene-SNPs


cobradata <- COBRAData(padj = results_padj)

cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = FDR)

basemethods(cobraperf)

colorscheme <- colors[basemethods(cobraperf)]

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colorscheme, incltruth = FALSE)


pdf(paste0(out_dir, "/venn_gene_snp_sign.pdf"))
plot_overlap(cobraplot, cex=c(1.2,1,0.7))
dev.off()


pdf(paste0(out_dir, "/upset_gene_snp_sign.pdf"))
plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 14, sets = basemethods(cobraperf), sets.bar.color = colors[basemethods(cobraperf)])
dev.off()



overlap <- cobraplot@overlap


summary$gene_snp_sign <- c(colSums(overlap), sum(rowSums(overlap == 1) == 2))



### Significant genes

# Keep minimal p-value per gene
library(matrixStats)

gene_ids <- factor(strsplit2(rownames(results_padj), ":")[, 1])

results_padj_split <- by(results_padj, gene_ids, function(x){
  
  colMins(as.matrix(x), na.rm = TRUE)
  
}, simplify = FALSE)


results_padj_gene <- data.frame(do.call(rbind, results_padj_split))
colnames(results_padj_gene) <- colnames(results_padj)

# Use iCOBRA
cobradata <- COBRAData(padj = results_padj_gene)

cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = FDR)

basemethods(cobraperf)

colorscheme <- colors[basemethods(cobraperf)]

cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colorscheme, incltruth = FALSE)


pdf(paste0(out_dir, "/venn_gene_sign.pdf"))
plot_overlap(cobraplot, cex=c(1.2,1,0.7))
dev.off()


pdf(paste0(out_dir, "/upset_gene_sign.pdf"))
plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 14, sets = basemethods(cobraperf), sets.bar.color = colors[basemethods(cobraperf)])
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
plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 14, sets = basemethods(cobraperf), sets.bar.color = colors[basemethods(cobraperf)])
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
plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 14, sets = basemethods(cobraperf), sets.bar.color = colors[basemethods(cobraperf)])
dev.off()


overlap <- cobraplot@overlap

summary$gene_all <- c(colSums(overlap), sum(rowSums(overlap == 1) == 2))




write.table(summary, paste0(out_dir, "/summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
















































