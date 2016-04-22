##############################################################################
## <<geuvadis_drimseq_0_3_3_comparison_plots.R>>

# BioC 3.2
# Created 2 Mar 2016 
# Modified 13 Apr 2016

# Plot the overlap versus top ranked sQTLs for drimseq and sqtlseeker

##############################################################################

Sys.time()

##############################################################################
# Libraries
##############################################################################

library(ggplot2)
library(plyr)


##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/geuvadis'
# population='CEU'
# comparison_out='drimseq_0_3_3_comparison_permutations'
# Overlaps_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotOverlaps.R'
# CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R'
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

print(rwd)
print(population)
print(Overlaps_function_path)
print(CAT_function_path)

##############################################################################

setwd(rwd)

comparison_out <- paste0(comparison_out, "/")

dir.create(comparison_out, showWarnings = FALSE, recursive = TRUE)


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


metadata <- rbind.fill(metadata)


metadata$method_name <- factor(metadata$method_name, levels = colors_df$methods)
metadata$method_name <- factor(metadata$method_name)


### gene_id column is used in calculateOverlaps as indicator of observations

results <- lapply(results, function(x){
  
  x$gene_id <- x$gene_snp
  
  return(x)
  
})



############################################################################
# Plot overlap versus top x DS genes
############################################################################


# source(Overlaps_function_path)
# 
# data_Overlaps <- list()
# 
# data_Overlaps[[1]] <- calculateOverlaps(results1 = results[["sqtlseeker"]], results2 = results[["drimseq"]], by = 100, FDR = FDR)
# 
# 
# 
# reference_method <- "sqtlseeker"
# 
# metadata_ov <- metadata[metadata$method_name == "drimseq", , drop = FALSE]
# ### Drop unnecessary levels
# metadata_ov$method_name <- factor(metadata_ov$method_name)
# 
# 
# ggp <- plotOverlaps(data_Overlaps, metadata = metadata_ov, plot_var = "method_name", facet_var = NULL, plot_colors = colors[levels(metadata_ov$method_name)], plotx = TRUE, reference_color = colors[reference_method])
# 
# 
# ggp <- ggp + 
#   # coord_cartesian(xlim = c(0, 500), ylim = c(0, 300)) +
#   ylab("Overlap with sqtlseeker") +
#   xlab("Number of top ranked sQTLs")
# 
# 
# pdf(paste0(comparison_out, "overlap_top_ranked_genes.pdf"), width = 7, height = 7)
# print(ggp)
# dev.off()
# 
# 
# ggp <- ggp + 
#   coord_cartesian(xlim = c(0, 40000), ylim = c(0, 40000)) +
#   ylab("Overlap with sqtlseeker") +
#   xlab("Number of top ranked sQTLs")
# 
# 
# pdf(paste0(comparison_out, "overlap_top_ranked_genes_zoom.pdf"), width = 7, height = 7)
# print(ggp)
# dev.off()
# 


############################################################################
# CAT plots - percentage overlap versus top x DS genes
############################################################################


source(CAT_function_path)

data_CAT <- list()

data_CAT[[1]] <- calculateCAT(results1 = results[["sqtlseeker"]], results2 = results[["drimseq"]], by = 100, FDR = FDR)

save(data_CAT, file = paste0(comparison_out, "data_CAT.Rdata"))


reference_method <- "sqtlseeker"

metadata_ov <- metadata[metadata$method_name == "drimseq", , drop = FALSE]
### Drop unnecessary levels
metadata_ov$method_name <- factor(metadata_ov$method_name)


ggp <- plotCAT(data_CAT, metadata = metadata_ov, plot_var = "method_name", facet_var = NULL, plot_colors = colors[levels(metadata_ov$method_name)], plotx = TRUE, reference_color = colors[reference_method])


ggp <- ggp + 
  coord_cartesian(ylim = c(0, 1)) +
  ylab("Percentage overlap with sqtlseeker") +
  xlab("Number of top ranked sQTLs")


pdf(paste0(comparison_out, "cat.pdf"), width = 7, height = 7)
print(ggp)
dev.off()


ggp <- ggp + 
  coord_cartesian(xlim = c(0, 10000), ylim = c(0, 1)) +
  ylab("Percentage overlap with sqtlseeker") +
  xlab("Number of top ranked sQTLs")


pdf(paste0(comparison_out, "cat_zoom.pdf"), width = 7, height = 7)
print(ggp)
dev.off()


ggp <- ggp + 
  coord_cartesian(xlim = c(0, 6000), ylim = c(0, 1)) +
  ylab("Percentage overlap with sqtlseeker") +
  xlab("Number of top ranked sQTLs")


pdf(paste0(comparison_out, "cat_zoom2.pdf"), width = 7, height = 7)
print(ggp)
dev.off()

save(ggp, file = paste0(comparison_out, "cat.Rdata"))

############################################################################
# CAT plots with CATplot from ffpe
############################################################################
# library(ffpe)

# vec1 <- results[["sqtlseeker"]]
# vec1 <- vec1[order(vec1$pvalue, decreasing = FALSE), "gene_snp"]
# 
# vec2 <- results[["drimseq"]]
# vec2 <- vec2[order(vec2$pvalue, decreasing = FALSE), "gene_snp"]
# 
# pdf(paste0(comparison_out, "cat_ffpe.pdf"), width = 7, height = 7)
# CATplot(vec1, vec2, maxrank = 10000, make.plot = TRUE)
# dev.off()












