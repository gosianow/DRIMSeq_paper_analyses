######################################################
## <<kim_drimseq_0_3_3_comparison_summary.R>>

# BioC 3.2
# Created 15 Jan 2015 
# Modified 14 Apr 2016

##############################################################################
Sys.time()
##############################################################################

library(plyr)
library(ggplot2)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/kim_adenocarcinoma'
# comparison_out='drimseq_0_3_3_comparison'
# keep_methods=c('dexseq','drimseq_genewise_grid_none','drimseq_genewise_grid_common','drimseq_genewise_grid_trended')
# text_size=18
# legend_size=16
# strip_size=16

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

dir.create(paste0(comparison_out, "figures/"), showWarnings = FALSE, recursive = TRUE)

out_dir <- paste0(comparison_out)

### colors

load(paste0(rwd, "/", comparison_out, "colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)


colors <- colors[keep_methods]
colors_df <- colors_df[colors_df$methods %in% keep_methods, , drop = FALSE]


##############################################################################



files <- list.files(comparison_out, pattern = "_summary.txt")
files


summary_list <- lapply(1:length(files), function(i){
  
  sm <- read.table(paste0(comparison_out, files[i]), header = TRUE, as.is = TRUE)
  
})


summary <- rbind.fill(summary_list)

write.table(summary, file = paste0(comparison_out, "summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




summary <- read.table(paste0(comparison_out, "summary.txt"), header = TRUE)

summary <- summary[summary$ds_method %in% colors_df$methods, , drop = FALSE]
summary$ds_method <- factor(summary$ds_method, levels = colors_df$methods)




ggp <- ggplot(data = summary, aes(x = count_method, y = counts_genes_ds, fill = ds_method)) +
  geom_bar(stat = "identity", position="dodge") +
  facet_wrap(~ model, nrow = 1) +
  ylab("Number of genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = text_size), axis.text.y = element_text(size = text_size), axis.title.x = element_blank(), axis.title.y = element_text(size = text_size, face = "bold"), legend.text = element_text(size = legend_size), legend.title = element_blank(), legend.position = "bottom", strip.text = element_text(size = strip_size)) +
  guides(fill = guide_legend(nrow = 1)) + 
  scale_fill_manual(values = colors[levels(summary$ds_method)])



pdf(paste0(comparison_out, "figures/summary_genes_ds_all.pdf"), 14, 7)
print(ggp)
dev.off()



ggp <- ggplot(data = summary, aes(x = count_method, y = counts_genes_all, fill = ds_method)) +
  geom_bar(stat = "identity", position="dodge") +
  facet_wrap(~ model, nrow = 1) +
  ylab("Number of genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = text_size), axis.text.y = element_text(size = text_size), axis.title.x = element_blank(), axis.title.y = element_text(size = text_size, face = "bold"), legend.text = element_text(size = legend_size), legend.title = element_blank(), legend.position = "bottom", strip.text = element_text(size = strip_size)) +
  guides(fill = guide_legend(nrow = 1)) + 
  scale_fill_manual(values = colors[levels(summary$ds_method)])



pdf(paste0(comparison_out, "figures/summary_genes_all.pdf"), 14, 7)
print(ggp)
dev.off()





summary_overlap <- summary[, c("model", "count_method", "ds_method", "counts_genes_ds_overlap")]
summary_overlap$set <- "overlap_with_dexseq"
colnames(summary_overlap)[4] <- "counts_genes_ds"


summary_unique <- summary[, c("model", "count_method", "ds_method", "counts_genes_ds")]
summary_unique$set <- "unique"
colnames(summary_unique)[4] <- "counts_genes_ds"


summary2 <- rbind(summary_overlap, summary_unique)

summary2$set <- factor(summary2$set, levels = c("overlap_with_dexseq", "unique"))


ggp <- ggplot(data = summary2, aes(x = count_method, y = counts_genes_ds, group = ds_method, fill = ds_method, alpha = set)) +
  geom_bar(stat = "identity", position="dodge") +
  facet_wrap(~ model, nrow = 1) +
  ylab("Number of genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = text_size), axis.text.y = element_text(size = text_size), axis.title.x = element_blank(), axis.title.y = element_text(size = text_size, face = "bold"), legend.text = element_text(size = legend_size), legend.title = element_blank(), legend.position = "bottom", strip.text = element_text(size = strip_size)) +
  guides(fill = guide_legend(nrow = 1)) + 
  scale_fill_manual(values = colors[levels(summary$ds_method)]) +
  scale_alpha_manual(values = c(1, 0.5))



pdf(paste0(comparison_out, "figures/summary_genes_ds_unique_and_overlap.pdf"), 14, 7)
print(ggp)
dev.off()





sessionInfo()










