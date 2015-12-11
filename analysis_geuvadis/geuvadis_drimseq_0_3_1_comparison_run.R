######################################################
## ----- geuvadis_drimseq_0_3_1_comparison_run
## <<geuvadis_drimseq_0_3_1_comparison_run.R>>

# BioC 3.1
# Created 18 Nov 2015 

##############################################################################

library(iCOBRA)
library(DRIMSeq)
library(ggplot2)
library(limma)
library(plyr)
library(Gviz)
library(GenomicFeatures)
library(tools)
library(rtracklayer)

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/Shared/data/seq/geuvadis'
population='CEU'


##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(population)


##############################################################################

setwd(rwd)

method_out_dir <- "drimseq_0_3_1_analysis/"

comparison_out <- "drimseq_0_3_1_comparison/"
dir.create(comparison_out, recursive = TRUE, showWarnings = FALSE)

out_dir <- comparison_out


##############################################################################

results <- list()


#####################################
### DRIMSeq results + adjust p-value
#####################################


res_list <- lapply(1:22, function(chr){
  # chr = 1
  
  res <- read.table(paste0(method_out_dir, population, "_chr",chr, "_results.txt"), header = TRUE, sep = "\t", as.is = TRUE)
  
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



#####################################
### sqtlseeker results
#####################################


res_tmp <- read.table(paste0("sqtlseeker_2_1_analysis/results/", population, "_results_all.txt"), header = TRUE, as.is = TRUE)
head(res_tmp)

colnames(res_tmp) <- c("gene_id", "snp_id", "F", "nb.groups", "md", "tr.first", "tr.second", "nb.perms", "pvalue")
res_tmp <- unique(res_tmp)


res_tmp$gene_snp <- paste0(res_tmp$gene_id, ":", res_tmp$snp_id)

res_tmp$adj_pvalue <- qvalue::qvalue(res_tmp$pvalue)$qvalues



results[["sqtlseeker"]] <- res_tmp

dim(results[["sqtlseeker"]])

table(duplicated(results[["sqtlseeker"]][, "gene_snp"]))


##########################################################################
### histograms of p-values
### scatterplot of p-values
##########################################################################


ggp <- DRIMSeq:::dm_plotPvalues(pvalues = results[["drimseq"]][, "pvalue"])

pdf(paste0(out_dir, "drimseq_hist_pvalues.pdf"))
print(ggp)
dev.off()


ggp <- DRIMSeq:::dm_plotPvalues(pvalues = results[["sqtlseeker"]][, "pvalue"])

pdf(paste0(out_dir, "sqtlseeker_hist_pvalues.pdf"))
print(ggp)
dev.off()



results_pval <- list()

pval_tmp <- results[["drimseq"]][, c("gene_snp", "pvalue")]
colnames(pval_tmp) <- c("gene_snp", "drimseq")
results_pval[["drimseq"]] <- pval_tmp

pval_tmp <- results[["sqtlseeker"]][, c("gene_snp", "pvalue")]
colnames(pval_tmp) <- c("gene_snp", "sqtlseeker")
results_pval[["sqtlseeker"]] <- pval_tmp


results_pval <- Reduce(function(...) merge(..., by = "gene_snp", all=TRUE, sort = FALSE), results_pval)
rownames(results_pval) <- results_pval$gene_snp
results_pval <- results_pval[, -1]


pdf(paste0(out_dir, "scatter_pvalues.pdf"))
smoothScatter(results_pval[, "drimseq"], results_pval[, "sqtlseeker"])
dev.off()


pdf(paste0(out_dir, "scatter_pvalues_range.pdf"))
smoothScatter(results_pval[, "drimseq"], results_pval[, "sqtlseeker"], xlim = c(0, 0.05))
dev.off()



# ggp <- ggplot(results_pval, aes(x = drimseq, y = sqtlseeker)) + 
#   stat_binhex(bins = 100) +
#   coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
# 
# pdf(paste0(out_dir, "binhex_pvalues.pdf"))
# print(ggp)
# dev.off()




##########################################################################
###  use iCOBRA for overlaps
##########################################################################

summary <- data.frame(method = c("drimseq", "sqtlseeker", "overlap"))


results_padj <- list()

padj_tmp <- results[["drimseq"]][, c("gene_snp", "adj_pvalue")]
colnames(padj_tmp) <- c("gene_snp", "drimseq")
results_padj[["drimseq"]] <- padj_tmp

padj_tmp <- results[["sqtlseeker"]][, c("gene_snp", "adj_pvalue")]
colnames(padj_tmp) <- c("gene_snp", "sqtlseeker")
results_padj[["sqtlseeker"]] <- padj_tmp


results_padj <- Reduce(function(...) merge(..., by = "gene_snp", all=TRUE, sort = FALSE), results_padj)
rownames(results_padj) <- results_padj$gene_snp
results_padj <- results_padj[, -1]



### COBRA

# All tested cases


cobradata <- COBRAData(padj = results_padj)

cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 1.1)

basemethods(cobraperf)


cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = c("darkorange1", "dodgerblue3"), incltruth = FALSE)

pdf(paste0(out_dir, "/venn_gene_snp_all.pdf"))
plot_overlap(cobraplot, cex=c(1.2,1,0.7))
dev.off()


overlap <- cobraplot@overlap

summary$gene_snp_all <- c(colSums(overlap), sum(rowSums(overlap == 1) == 2))

overlap <- cobraplot@overlap[!duplicated(strsplit2(rownames(cobraplot@overlap), ":")[, 1]), ]

summary$gene_all <- c(colSums(overlap), sum(rowSums(overlap == 1) == 2))


pdf(paste0(out_dir, "/venn_gene_all.pdf"))

vennDiagram(overlap, circle.col = c("darkorange1", "dodgerblue3"))

dev.off()



# Significant cases


cobradata <- COBRAData(padj = results_padj)

cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 0.05)

basemethods(cobraperf)


cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = c("darkorange1", "dodgerblue3"), incltruth = FALSE)

pdf(paste0(out_dir, "/venn_gene_snp_sign.pdf"))
plot_overlap(cobraplot, cex=c(1.2,1,0.7))
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

vennDiagram(overlap, circle.col = c("darkorange1", "dodgerblue3"))

dev.off()




write.table(summary, paste0(out_dir, "/summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




##########################################################################
### Load list of validated genes
##########################################################################


valid <- read.table("data/validation/positive_controls_sqtlseeker_paper.txt", sep = "\t", header = TRUE, as.is = TRUE)

valid$gene_snp <- paste0(valid$gene_id, ":", valid$snp_id)
valid$chr <- strsplit2(valid$snp_id, "_")[, 2]

valid$snp_position <- as.numeric(strsplit2(valid$snp_id, "_")[, 3])



##########################################################################
### Proportion plots of validated sQTLs
##########################################################################

### drimseq plots

valid$gene_snp %in% results[["drimseq"]]$gene_snp

results[["drimseq"]][which(results[["drimseq"]]$gene_snp %in% valid$gene_snp), ]


for(i in 1:nrow(valid)){
  # i = 1
  
  load(paste0(method_out_dir, population, "_chr",valid$chr[i], "_d.Rdata"))
  
  plot_main <- paste0(valid$sqtlseeker_gene_id[i], " - ", valid$sqtlseeker_snp_id[i], "\n FDR = ", sprintf( "%.02e",results[["drimseq"]][which(results[["drimseq"]]$gene_snp == valid$gene_snp[i]), "adj_pvalue"]))
  
  plotFit(d, gene_id = valid$gene_id[i], snp_id = valid$snp_id[i], plot_type = "boxplot1", order = FALSE, plot_full = TRUE, plot_main = plot_main, out_dir = paste0(out_dir, "/drimseq1_"))
  
  plotFit(d, gene_id = valid$gene_id[i], snp_id = valid$snp_id[i], plot_type = "boxplot2", order = FALSE, plot_full = TRUE, plot_main = plot_main, out_dir = paste0(out_dir, "/drimseq2_"))
  
  plotFit(d, gene_id = valid$gene_id[i], snp_id = valid$snp_id[i], plot_type = "boxplot1", order = TRUE, plot_full = TRUE, plot_main = plot_main, out_dir = paste0(out_dir, "/drimseq1o_"))
  
  plotFit(d, gene_id = valid$gene_id[i], snp_id = valid$snp_id[i], plot_type = "boxplot2", order = TRUE, plot_full = TRUE, plot_main = plot_main, out_dir = paste0(out_dir, "/drimseq2o_"))
  
  
}



### sqtlseeker plots

valid$gene_snp %in% results[["sqtlseeker"]]$gene_snp

results[["sqtlseeker"]][which(results[["sqtlseeker"]]$gene_snp %in% valid$gene_snp), ]


sqtlseeker_counts <- read.table("sqtlseeker_2_1_analysis/data/trExpCount_CEU_sqtlseeker_ratios.tsv", header = TRUE, sep = "\t", as.is = TRUE)


for(i in 1:nrow(valid)){
  # i = 1
  
  load(paste0(method_out_dir, population, "_chr",valid$chr[i], "_d.Rdata"))
  
  counts <- as.matrix(sqtlseeker_counts[sqtlseeker_counts$geneId == valid$gene_id[i], -c(1, 2)])
  rownames(counts) <- sqtlseeker_counts[sqtlseeker_counts$geneId == valid$gene_id[i], "trId"]
  
  block <- d@blocks[[valid$gene_id[i]]]
  block_id <- block[which(block[, "snp_id"] == valid$snp_id[i]), "block_id"]
  
  group <- factor(d@genotypes[[valid$gene_id[i]]][block_id, colnames(counts)])
  
  plot_main <- paste0(valid$sqtlseeker_gene_id[i], " - ", valid$sqtlseeker_snp_id[i], "\n FDR = ", sprintf( "%.02e",results[["sqtlseeker"]][which(results[["sqtlseeker"]]$gene_snp == valid$gene_snp[i]), "adj_pvalue"]))
  
  
  ggp <- DRIMSeq:::dm_plotProportions(counts, group, pi_full = NULL, pi_null = NULL, main = plot_main, plot_type = "boxplot1", order = FALSE)
  
  
  pdf(paste0(out_dir, "/sqtlseeker1_",  gsub(pattern = "\\.", replacement = "_" , paste0( valid$gene_id[i], "_", valid$snp_id[i])), ".pdf"), width = 12, height = 7)
  print(ggp)
  dev.off()
  
  
  ggp <- DRIMSeq:::dm_plotProportions(counts, group, pi_full = NULL, pi_null = NULL, main = plot_main, plot_type = "boxplot2", order = FALSE)
  
  
  pdf(paste0(out_dir, "/sqtlseeker2_",  gsub(pattern = "\\.", replacement = "_" , paste0( valid$gene_id[i], "_", valid$snp_id[i])), ".pdf"), width = 12, height = 7)
  print(ggp)
  dev.off()
  
  
}




##########################################################################
### Gene structure plot for validated snps
##########################################################################

path_gtf <- "geuvadis_annotation/gencode.v12.annotation.gtf"

gat <- GenomeAxisTrack()

gtf <- import(path_gtf)

txdb <- makeTxDbFromGFF(path_gtf, format="gtf")
genes <- genes(txdb)


txTr <- GeneRegionTrack(txdb, name = "transcripts",  transcriptAnnotation = "transcript", just.group = "above", cex.group = 1)


for(j in 1:nrow(valid)){
  # j = 1
  
  g <- valid$gene_id[j]
  chromosome <- as.character(seqnames(genes[g, ]))
  astart <- start(genes[g, ]) - 5000
  aend <- end(genes[g, ]) + 5000
  
  
  # gtf_sub <- gtf[mcols(gtf)$gene_id == g, ]
  # 
  # txTr <- GeneRegionTrack(gtf_sub, name = "transcripts",  transcriptAnnotation = "transcript", just.group = "above", cex.group = 1)
  
  
  ht <- HighlightTrack(trackList = list(gat, txTr), start = c(valid$snp_position[j]), width = 1, chromosome = chromosome)
  
  
  
  pdf(paste0(out_dir, "/gviz_", valid$sqtlseeker_gene_id[j], ".pdf"), 7, 4)
  
  plotTracks(c(ht), from = astart, to = aend, chromosome = chromosome)
  
  dev.off()
  
  
  
  pdf(paste0(out_dir, "/gviz_", valid$sqtlseeker_gene_id[j], "_zoom.pdf"), 7, 4)
  
  plotTracks(c(ht), from = valid$snp_position[j] - 5000, to = valid$snp_position[j] + 5000, chromosome = chromosome)
  
  dev.off()
  
}

















































