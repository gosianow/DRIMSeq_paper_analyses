library(ggplot2)


methods <- factor(c("dexseq","drimseq_common", "drimseq_genewise_grid_common", "drimseq_genewise_grid_none", "drimseq_genewise_constrOptim", "drimseq_genewise_optim", "drimseq_genewise_optimize"), levels = c("dexseq","drimseq_common", "drimseq_genewise_grid_common", "drimseq_genewise_grid_none", "drimseq_genewise_constrOptim", "drimseq_genewise_optim", "drimseq_genewise_optimize"))


colors_df <- data.frame(methods = methods, colors = c("#4065B1", "#7FB972", "#D92120", "#E68E34", "#C71585", "#781C81", "#B17BA6"))
colors_df$colors <- as.character(colors_df$colors)

ggp <- ggplot(colors_df, aes(x = methods, y=rep(1, length(colors)), fill = methods)) + geom_bar(stat="identity") + scale_fill_manual(values = colors_df$colors) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


pdf(paste0("/home/gosia/R/drimseq_paper/colors.pdf"), 10, 5)
print(ggp)
dev.off()


colors <- colors_df$colors
names(colors) <- colors_df$methods


write.table(colors_df, paste0("/home/gosia/R/drimseq_paper/colors.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(colors, colors_df, file = paste0("/home/gosia/R/drimseq_paper/colors.Rdata"))




write.table(colors_df, paste0("/home/gosia/multinomial_project/simulations_sim5/colors.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(colors, colors_df, file = paste0("/home/gosia/multinomial_project/simulations_sim5/colors.Rdata"))



write.table(colors_df, paste0("/home/Shared/data/seq/kim_adenocarcinoma/drimseq_0_3_1_comparison/colors.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(colors, colors_df, file = paste0("/home/Shared/data/seq/kim_adenocarcinoma/drimseq_0_3_1_comparison/colors.Rdata"))



write.table(colors_df, paste0("/home/Shared/data/seq/brooks_pasilla/drimseq_0_3_1_comparison/colors.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

save(colors, colors_df, file = paste0("/home/Shared/data/seq/brooks_pasilla/drimseq_0_3_1_comparison/colors.Rdata"))








