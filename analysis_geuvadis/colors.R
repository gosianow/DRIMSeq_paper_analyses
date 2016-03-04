##############################################################################
## <<colors.R>>

# BioC 3.2
# Created 1 Mar 2016 

##############################################################################

Sys.time()

##############################################################################

library(ggplot2)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/geuvadis'
# out_dir='/home/Shared/data/seq/geuvadis/drimseq_0_3_3_comparison'

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
print(out_dir)

##############################################################################

setwd(rwd)

dir.create(out_dir, recursive = TRUE)

##############################################################################


methods <- factor(c("sqtlseeker","drimseq"), levels = c("sqtlseeker","drimseq"))


colors_df <- data.frame(methods = methods, colors = c("forestgreen", "orange2"), colors_hex = c("#228B22", "#EE9A00"))
colors_df$colors <- as.character(colors_df$colors)
colors_df$colors_hex <- as.character(colors_df$colors_hex)


ggp <- ggplot(colors_df, aes(x = methods, y=rep(1, length(colors)), fill = methods)) + geom_bar(stat="identity") + scale_fill_manual(values = colors_df$colors) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


pdf(paste0(out_dir, "/colors.pdf"), 10, 5)
print(ggp)
dev.off()


colors <- colors_df$colors
names(colors) <- colors_df$methods



write.table(colors_df, paste0(out_dir, "/colors.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(colors, colors_df, file = paste0(out_dir, "/colors.Rdata"))








