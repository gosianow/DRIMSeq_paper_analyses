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


res_tmp <- lapply(1:22, function(chr){
  # chr = 1
  
  res <- read.table(paste0(method_out_dir, population, "_chr",chr, "_results.txt"), header = TRUE, sep = "\t", as.is = TRUE)
  
  return(res)
  
})


res <- rbind.fill(res_tmp)

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


res_tmp <- read.table(paste0("sqtlseeker_2_0_analysis/Results/", population, "_results_all.txt"), header = TRUE, as.is = TRUE)
colnames(res_tmp) <- c("gene_id", "snp_id", "F", "nb.groups", "md", "tr.first", "tr.second", "nb.perms", "pvalue")
res_tmp <- unique(res_tmp)


res_tmp$gene_snp <- paste0(res_tmp$gene_id, ":", res_tmp$snp_id)

res_tmp$qvalue <- qvalue::qvalue(res_tmp$pvalue)$qvalues

res_tmp$adj_pvalue <- p.adjust(res_tmp$pvalue, method = "BH")

table(res_tmp$qvalue < 0.05)
# table(res_tmp$qvalue < 0.05 & res_tmp$md >= 0.01)

table(res_tmp$adj_pvalue < 0.05)


results[["sqtlseeker"]] <- res_tmp

dim(results[["sqtlseeker"]])

table(duplicated(results[["sqtlseeker"]][, "gene_snp"]))


##########################################################################
### histograms of p-values
##########################################################################



DRIMSeq:::dm_plotTable(pvalues = results[["drimseq"]][, "pvalue"], name = "pvalue", out_dir = paste0(out_dir, "drimseq_"))

DRIMSeq:::dm_plotTable(pvalues = results[["sqtlseeker"]][, "pvalue"], name = "pvalue", out_dir = paste0(out_dir, "sqtlseeker_"))



##########################################################################
###  use iCOBRA
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
































