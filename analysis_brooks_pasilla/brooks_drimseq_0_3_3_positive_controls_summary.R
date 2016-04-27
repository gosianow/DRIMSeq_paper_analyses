######################################################
## ----- brooks_drimseq_0_3_3_positive_controls_summary
## <<brooks_drimseq_0_3_3_positive_controls_summary.R>>

# BioC 3.2
# Created 15 Jan 2015 

##############################################################################

Sys.time()

##############################################################################

library(plyr)
library(ggplot2)
library(rtracklayer)
library(reshape2)
library(Gviz)
library(GenomicFeatures)
library(tools)
library(limma)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/brooks_pasilla'
# method_out='drimseq_0_3_3'
# comparison_out='drimseq_0_3_3_positive_controls'
# keep_methods=c('dexseq','drimseq_genewise_grid_none','drimseq_genewise_grid_common','drimseq_genewise_grid_trended')

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

dir.create(paste0(comparison_out, "figures/"), showWarnings = FALSE, recursive = TRUE)


##############################################################################
# merge positive controls summary into one table
##############################################################################


files <- list.files(comparison_out, pattern = "_validation_summary.txt")
files


summary_list <- lapply(1:length(files), function(i){
  
  sm <- read.table(paste0(comparison_out, files[i]), header = TRUE, as.is = TRUE)
  
})


summary <- rbind.fill(summary_list)

summary <- summary[, c("model", "count_method", keep_methods)]


write.table(summary, file = paste0(comparison_out, "validation_summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



summarym <- melt(summary, id.vars = c("model", "count_method"))
# summarym <- melt(summary, id.vars = c("model", "count_method"), variable.name = "ds_method", value.name = "counts")

summarym$variable <- factor(summarym$variable, levels = rev(keep_methods))


ggp <- ggplot(summarym, aes(x = count_method, y = variable, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = value), color = "black", size = 7) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  xlab("") + 
  ylab("") + 
  theme_bw() +
  theme(panel.background = element_rect(fill = NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14), strip.text = element_text(size = 14)) +
  scale_fill_gradient(low = "grey90", high = "grey50", na.value = "grey90") +
  facet_wrap(~ model)


pdf(paste0(comparison_out, "figures/validation_summary.pdf"), 15, 5)
print(ggp)
dev.off()



##############################################################################
# merge positive controls into tables per model
##############################################################################


files <- list.files(comparison_out, pattern = "_validation.txt")
files


summary_list <- lapply(1:length(files), function(i){
  
  sm <- read.table(paste0(comparison_out, files[i]), header = TRUE, as.is = TRUE)
  
})


summary <- rbind.fill(summary_list)

summary <- summary[, c("brooks_gene_id", "gene_id", "model", "count_method", keep_methods)]
summary$model <- factor(summary$model)

models <- levels(summary$model)




for(i in 1:nlevels(summary$model)){
  # i = 1
  
  # summarym <- melt(summary[summary$model == models[i], ], id.vars = c("brooks_gene_id", "gene_id", "model", "count_method"), variable.name = "ds_method", value.name = "counts")
  summarym <- melt(summary[summary$model == models[i], ], id.vars = c("brooks_gene_id", "gene_id", "model", "count_method"))
  
  summarym$variable <- factor(summarym$variable, levels = keep_methods)
  
  summarym$status <- factor(summarym$value < 0.05)
  summarym$full_gene_id <- paste0(summarym$brooks_gene_id, "_", summarym$gene_id)
  
  ggp <- ggplot(summarym, aes(x = variable, y = full_gene_id, fill = status)) + 
    geom_tile() + 
    geom_text(aes(label = sprintf( "%.02e", value)), color = "black", size = 4) + 
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0)) + 
    xlab("") + 
    ylab("") + 
    theme_bw() +
    theme(panel.background = element_rect(fill = NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14), strip.text = element_text(size = 14)) +
    scale_fill_manual(values = c("grey80", "grey50"), na.value = "grey90") +
    facet_wrap(~ count_method, nrow = 1)
  
  
  pdf(paste0(comparison_out, "figures/", models[i], "_validation_summary.pdf"), 16, 10)
  print(ggp)
  dev.off()
  
  
}




sessionInfo()

