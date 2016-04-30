##############################################################################
# <<geuvadis_drimseq_0_3_3_comparison_run.R>>

# BioC 3.2
# Created 29 Feb 2016 

# Merge results and recalculate adjusted p-values for drimseq and sqtlseeker
# Create a table with results for validated sQTLs
# Create a summary table with number of detected validated sQTLs
# Plot proportions for validated sQTLs

##############################################################################

Sys.time()

##############################################################################

library(DRIMSeq)
library(reshape2)
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

# rwd='/home/Shared/data/seq/geuvadis'
# population='CEU'
# valid_path='data/validation/glimmps/glimmps_valid_pcr.txt'
# plot_proportions=TRUE
# plot_tables=TRUE
# method_out='drimseq_0_3_3_analysis_permutations_all_genes'
# comparison_out='drimseq_0_3_3_comparison_permutations_all_genes'
# positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes'
# sqtlseeker_results='sqtlseeker_2_1_analysis'
# FDR=0.05


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

positive_controls_out <- paste0(positive_controls_out, "/", basename(file_path_sans_ext(valid_path)), "/")

dir.create(positive_controls_out, recursive = TRUE, showWarnings = FALSE)

out_dir <- positive_controls_out

dir.create(paste0(out_dir, "figures/"), recursive = TRUE, showWarnings = FALSE)

comparison_out <- paste0(comparison_out, "/")
method_out <- paste0(method_out, "/")



##############################################################################
### colors
##############################################################################

load(paste0(rwd, "/", comparison_out, "colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)

##############################################################################
# validated genes
##############################################################################


valid <- read.table(valid_path, header = TRUE, sep = "\t", as.is = TRUE) 


##############################################################################
# merge results 
##############################################################################


results <- list()

#####################################
### sqtlseeker results
#####################################


res <- read.table(paste0(comparison_out, "results_sqtlseeker.txt"), header = TRUE, as.is = TRUE)
head(res)


results[["sqtlseeker"]] <- res

#####################################
### DRIMSeq results
#####################################


res <- read.table(paste0(comparison_out, "results_drimseq.txt"), header = TRUE, as.is = TRUE)
head(res)

results[["drimseq"]] <- res


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


############################################################################
# Check which validated sqtls are significant 
############################################################################


results_padj_valid <- data.frame(matrix(NA, nrow = nrow(valid), ncol = ncol(results_padj)))
rownames(results_padj_valid) <- valid$gene_snp
colnames(results_padj_valid) <- colnames(results_padj)

mm <- match(rownames(results_padj_valid), rownames(results_padj))

results_padj_valid[,] <- results_padj[mm, ]

head(results_padj_valid)


out <- data.frame(valid[, c("gene_id", "gene_name", "snp_id", "snp_name", "gene_snp")], results_padj_valid, stringsAsFactors = FALSE)

write.table(out, file = paste0(positive_controls_out, "validation.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




number_sign_valid <- matrix(colSums(results_padj_valid < FDR, na.rm = TRUE), nrow = 1)
colnames(number_sign_valid) <- colnames(results_padj_valid)


out_summary <- data.frame(number_sign_valid, validated_tested = sum(!is.na(mm)), validated = nrow(results_padj_valid))
out_summary

write.table(out_summary, file = paste0(positive_controls_out, "validation_summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



##############################################################################
# Plot tables with validated sQTLs
##############################################################################

keep_methods <- colors_df$methods

if(plot_tables){
  
  summary <- read.table(paste0(positive_controls_out, "validation.txt"), header = TRUE)
  
  summarym <- melt(summary, id.vars = c("gene_id", "gene_name", "snp_id", "snp_name", "gene_snp"))
  
  summarym$variable <- factor(summarym$variable, levels = keep_methods)
  summarym$status <- factor(summarym$value < FDR)
  
  summarym$gene_snp_name <- factor(paste0(summarym$gene_name, " : ", summarym$snp_name), levels = paste0(summary$gene_name, " : ", summary$snp_name)[nrow(summary):1])
  
  
  ggp <- ggplot(summarym, aes(x = variable, y = gene_snp_name, fill = status)) + 
    geom_tile() + 
    geom_text(aes(label = sprintf( "%.02e", value)), color = "black", size = 4) + 
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0)) + 
    xlab("") + 
    ylab("") + 
    theme_bw() +
    theme(panel.background = element_rect(fill = NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(size = 14, angle = 0, vjust = 0, hjust = 0.5), axis.text.y = element_text(size = 14), strip.text = element_text(size = 14), legend.position = "none") +
    scale_fill_manual(values = c("grey80", "grey50"), na.value = "grey90")
  
  
  pdf(paste0(positive_controls_out, "validation.pdf"), 7, 7)
  print(ggp)
  dev.off()
  
  
}


##########################################################################
### Proportion plots of validated sQTLs
##########################################################################

# library(devtools)
# load_all("/home/gosia/R/package_devel/DRIMSeq")

if(plot_proportions){
  
  
  ### drimseq plots
  
  for(i in 1:length(valid$gene_snp)){
    # i = 8
    
    load(paste0(method_out, population, "_chr",valid$chr[i], "_d.Rdata"))
    
    if(!valid$gene_snp %in% d@results$gene_snp)
      next
    
    mean_expression <- round(d@mean_expression[valid$gene_id[i]], 1)
    adj_pvalue <- results[["drimseq"]][which(results[["drimseq"]]$gene_snp == valid$gene_snp[i]), "adj_pvalue"]
    
    plot_main <- paste0(valid$gene_name[i], " - ", valid$snp_name[i], "\n mean gene expression = ", mean_expression, " / FDR = ", sprintf( "%.02e", adj_pvalue))
    

    ggp <- plotFit(d, gene_id = valid$gene_id[i], snp_id = valid$snp_id[i], plot_type = "boxplot1", order = FALSE)
    
    ggp <- ggp + ggtitle(plot_main)
    
    pdf(paste0(out_dir, "figures/drimseq_boxplot1_", i, "_", valid$gene_name[i], "_", valid$snp_name[i], ".pdf"), width = 14, height = 7)
    print(ggp)
    dev.off()
    
    
    ggp <- plotFit(d, gene_id = valid$gene_id[i], snp_id = valid$snp_id[i], plot_type = "boxplot1", order = TRUE)
    
    ggp <- ggp + ggtitle(plot_main)
    
    pdf(paste0(out_dir, "figures/drimseq_boxplot1order_", i, "_", valid$gene_name[i], "_", valid$snp_name[i], ".pdf"), width = 14, height = 7)
    print(ggp)
    dev.off()
    
    
    ggp <- plotFit(d, gene_id = valid$gene_id[i], snp_id = valid$snp_id[i], plot_type = "boxplot3", order = TRUE)
    
    ggp <- ggp + ggtitle(plot_main)
    
    pdf(paste0(out_dir, "figures/drimseq_boxplot3order_", i, "_", valid$gene_name[i], "_", valid$snp_name[i], ".pdf"), width = 14, height = 7)
    print(ggp)
    dev.off()
    
    
    
  }
  
  
}




















