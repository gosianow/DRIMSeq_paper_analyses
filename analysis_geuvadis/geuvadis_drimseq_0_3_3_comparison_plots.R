##############################################################################
## <<geuvadis_drimseq_0_3_3_comparison_plots.R>>

# BioC 3.2
# Created 2 Mar 2016 

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
# Overlaps_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotOverlaps.R'
# CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R'

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

method_out <- "drimseq_0_3_3_analysis/"

comparison_out <- "drimseq_0_3_3_comparison/"
dir.create(comparison_out, showWarnings = FALSE, recursive = TRUE)


##############################################################################
# colors
##############################################################################


load(paste0(rwd, "/", "drimseq_0_3_3_comparison", "/colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)



#######################################################
# merge results 
#######################################################


##############################################################################
# merge results 
##############################################################################

results <- list()
metadata <- list()


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
metadata[["sqtlseeker"]] <- data.frame(method_name = "sqtlseeker", stringsAsFactors = FALSE)


# dim(results[["sqtlseeker"]])
# 
# table(duplicated(results[["sqtlseeker"]][, "gene_snp"]))



#####################################
### DRIMSeq results + adjust p-value
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
# data_Overlaps[[1]] <- calculateOverlaps(results1 = results[["sqtlseeker"]], results2 = results[["drimseq"]], by = 100)
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

data_CAT[[1]] <- calculateCAT(results1 = results[["sqtlseeker"]], results2 = results[["drimseq"]], by = 100)



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
  coord_cartesian(xlim = c(0, 40000), ylim = c(0, 1)) +
  ylab("Percentage overlap with sqtlseeker") +
  xlab("Number of top ranked sQTLs")


pdf(paste0(comparison_out, "cat_zoom.pdf"), width = 7, height = 7)
print(ggp)
dev.off()



















