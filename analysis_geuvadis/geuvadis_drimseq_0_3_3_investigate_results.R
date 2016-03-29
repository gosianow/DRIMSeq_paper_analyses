
# R32

# library(devtools)
# load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")

##############################################################################
# GEUVADIS 
##############################################################################

library(DRIMSeq)
library(ggplot2)
library(limma)
library(GenomicRanges)
library(plyr)

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/Shared/data/seq/geuvadis'
workers=10
population='CEU'
chr='1'


setwd(rwd)

out_dir <- "drimseq_0_3_3_analysis/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_name <- paste0(out_dir, population, "_chr",chr, "_")

load(paste0(out_name, "d.Rdata"))


### New outputdirectory
out_dir <- "drimseq_0_3_3_analysis/investigate_results/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_name <- paste0(out_dir, population, "_chr", chr, "_")


#####################################
### sqtlseeker results
#####################################


resseeker <- read.table(paste0("sqtlseeker_2_1_analysis/results/", population, "_results_all.txt"), header = TRUE, as.is = TRUE)
head(resseeker)

colnames(resseeker) <- c("gene_id", "snp_id", "F", "nb.groups", "md", "tr.first", "tr.second", "nb.perms", "pvalue")
resseeker <- unique(resseeker)


resseeker$gene_snp <- paste0(resseeker$gene_id, ":", resseeker$snp_id)

resseeker$adj_pvalue <- qvalue::qvalue(resseeker$pvalue)$qvalues



#####################################
### DRIMSeq results + adjust p-value
#####################################


res_list <- lapply(1:22, function(chr){
  # chr = 1
  
  res_tmp <- read.table(paste0(method_out, population, "_chr",chr, "_results.txt"), header = TRUE, sep = "\t", as.is = TRUE)
  
  return(res_tmp)
  
})


res <- rbind.fill(res_list)

res <- res[!is.na(res$pvalue), ]

res$gene_block <- paste0(res$gene_id, ":", res$block_id)
res$gene_snp <- paste0(res$gene_id, ":", res$snp_id)


res_uniq <- res[!duplicated(res$gene_block), ]
res_uniq$adj_pvalue <- p.adjust(res_uniq$pvalue, method = "BH")
mm <- match(res$gene_block, res_uniq$gene_block)
res$adj_pvalue <- res_uniq$adj_pvalue[mm]

resdrim <- res


####################################
# Plots of dispersion and df vs p-values
####################################


res <- results(d)

gene_ids <- names(d@genewise_dispersion)


genewise_disp <- lapply(1:length(d@genewise_dispersion), function(i){
  
  data.frame(gene_id = gene_ids[i], block_id = names(d@genewise_dispersion[[i]]), genewise_dispersion = d@genewise_dispersion[[i]], stringsAsFactors = FALSE)
  
})


genewise_disp <- rbind.fill(genewise_disp)


resm <- merge(res, genewise_disp, by = c("gene_id", "block_id"), sort = FALSE)

### Keep only unique snps - blocks
resm <- resm[!duplicated(resm[, c("gene_id", "block_id")]), ]



ggp <- ggplot(data = resm, aes(x = pvalue, y = log10(genewise_dispersion))) + 
  geom_hex(binwidth = c(0.01, 0.05)) +
  theme_bw() +
  ylab("gamma_+") +
  xlab("P-values") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6)))


pdf(paste0(out_name, "dispersion_vs_pvalues.pdf"))
print(ggp)
dev.off()


max_pval <- 0.05

ggp <- ggplot(data = resm[resm$pvalue < max_pval, ], aes(x = pvalue, y = log10(genewise_dispersion))) + 
  geom_hex(binwidth = c(0.0005, 0.05)) +
  theme_bw() +
  ylab("gamma_+") +
  xlab("P-values") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6)))


pdf(paste0(out_name, "dispersion_vs_pvalues_zoom.pdf"))
print(ggp)
dev.off()



ggp <- ggplot(data = resm, aes(x = log10(pvalue), y = log10(genewise_dispersion))) + 
  geom_hex(binwidth = c(0.01, 0.05)) +
  theme_bw() +
  ylab("gamma_+") +
  xlab("P-values") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6)))


pdf(paste0(out_name, "dispersion_vs_pvalues_log.pdf"))
print(ggp)
dev.off()



ggp <- ggplot(data = resm, aes(x = pvalue, y = df)) + 
  geom_hex(binwidth = c(0.01, 1)) +
  theme_bw() +
  ylab("Degrees of freedom") +
  xlab("P-values") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6)))


pdf(paste0(out_name, "df_vs_pvalues.pdf"))
print(ggp)
dev.off()


max_pval <- 0.05

ggp <- ggplot(data = resm[resm$pvalue < max_pval, ], aes(x = pvalue, y = df)) + 
  geom_hex(binwidth = c(0.0005, 1)) +
  theme_bw() +
  ylab("Degrees of freedom") +
  xlab("P-values") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6))) 


pdf(paste0(out_name, "df_vs_pvalues_zoom.pdf"))
print(ggp)
dev.off()



####################################
# Fraction of sign blocks for a gene
####################################

res <- results(d)

### Keep only unique snps - blocks
res <- res[!duplicated(res[, c("gene_id", "block_id")]), ]


res$gene_id <- factor(res$gene_id, levels = unique(res$gene_id))


frac_sign_blocks <- by(res, res$gene_id, function(x){
  
  mean(x$adj_pvalue < 0.05, na.rm = TRUE)
  
}, simplify = FALSE)


frac_sign_blocks <- data.frame(frac_sign_blocks = unlist(frac_sign_blocks))

frac_sign_blocks <- frac_sign_blocks[frac_sign_blocks$frac_sign_blocks > 0, , drop = FALSE]


ggp <- ggplot(frac_sign_blocks, aes(x = frac_sign_blocks)) +
  geom_histogram(binwidth = 0.01, fill = "grey") +
  theme_bw()


pdf(paste0(out_name, "frac_sign_blocks.pdf"))
print(ggp)
dev.off()



count_sign_blocks <- by(res, res$gene_id, function(x){
  
  sum(x$adj_pvalue < 0.05, na.rm = TRUE)
  
}, simplify = FALSE)


count_sign_blocks <- data.frame(count_sign_blocks = unlist(count_sign_blocks), gene_id = levels(res$gene_id), stringsAsFactors = FALSE)

count_sign_blocks <- count_sign_blocks[count_sign_blocks$count_sign_blocks > 0, , drop = FALSE]


ggp <- ggplot(count_sign_blocks, aes(x = count_sign_blocks)) +
  geom_histogram(binwidth = 1, fill = "grey") +
  theme_bw()


pdf(paste0(out_name, "count_sign_blocks.pdf"))
print(ggp)
dev.off()




### For chr1 there is a gene with 144 significant blocks (ENSG00000143799.7)

indx_max <- which.max(count_sign_blocks$count_sign_blocks)

gene_max <- count_sign_blocks[indx_max, "gene_id"]


genes <- res[res$gene_id == gene_max & res$adj_pvalue < 0.05, ]

genes <- genes[order(genes$pvalue, decreasing = FALSE), ]


genotypes <- d@genotypes[[gene_max]][genes$block_id, ]


### Plot a heat map with genotypes
library(pheatmap)
library(RColorBrewer)
library(colorRamps)


ph <- pheatmap(genotypes, cluster_cols = FALSE, fontsize = 6, filename = paste0(out_name, "genotypes_", gene_max ,".pdf"), width = 10, height = 14)


block_idph$tree_row$labels[ph$tree_row$order]



ph <- pheatmap(genotypes, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 6, filename = paste0(out_name, "genotypes_", gene_max ,"_sortsign.pdf"), width = 10, height = 14)





for(i in 1:10){
  
  plotFit(d, gene_id = genes[i, "gene_id"], snp_id = genes[i, "snp_id"], out_dir = paste0(out_name, i, "_"))
  
}






























