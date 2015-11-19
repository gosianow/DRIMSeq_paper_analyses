######################################################
## ----- kim_summary
## <<kim_summary.R>>

# BioC 3.1
# Created 18 Nov 2015 

##############################################################################

library(plyr)
library(ggplot2)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/kim_adenocarcinoma/'


##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)



##############################################################################

setwd(rwd)

comparison_out <- "drimseq_0_3_1_comparison/"


out_dir <- paste0(comparison_out)

### Colors

colors_tmp <- read.table(paste0(comparison_out, "colors.txt"), header = TRUE, as.is = TRUE)
colors <- colors_tmp$colors
names(colors) <- colors_tmp$methods


##############################################################################



files <- list.files(comparison_out, pattern = "summary.txt")
files


summary_list <- lapply(1:length(files), function(i){
  
  sm <- read.table(paste0(comparison_out, files[i]), header = TRUE, as.is = TRUE)
  
})


summary <- rbind.fill(summary_list)

write.table(summary, file = paste0(comparison_out, "summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




summary <- read.table(paste0(comparison_out, "summary.txt"), header = TRUE)



ggp <- ggplot(data = summary, aes(x = count_method, y = counts_genes_ds, fill = ds_method)) +
  geom_bar(stat = "identity", position="dodge") +
  facet_wrap(~ model, nrow = 1) +
  ylab("Number of genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16), axis.text.y = element_text(size = 16), axis.title.x = element_blank(), axis.title.y = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16), legend.title = element_blank(), legend.position = "bottom", strip.text = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 2)) + 
  scale_fill_manual(values = colors[levels(summary$ds_method)])



pdf(paste0(comparison_out, "summary_genes_ds.pdf"), 17, 7)
print(ggp)
dev.off()





ggp <- ggplot(data = summary, aes(x = count_method, y = counts_genes_all, fill = ds_method)) +
  geom_bar(stat = "identity", position="dodge") +
  facet_wrap(~ model, nrow = 1) +
  ylab("Number of genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16), axis.text.y = element_text(size = 16), axis.title.x = element_blank(), axis.title.y = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16), legend.title = element_blank(), legend.position = "bottom", strip.text = element_text(size = 14)) +
  guides(fill = guide_legend(nrow = 2)) + 
  scale_fill_manual(values = colors[levels(summary$ds_method)])



pdf(paste0(comparison_out, "summary_genes_all.pdf"), 14, 7)
print(ggp)
dev.off()




















