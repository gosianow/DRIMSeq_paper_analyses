######################################################
## ----- sim5_drimseq_comparison
## <<sim5_drimseq_comparison.R>>

# BioC 3.2
# Created 16 Nov 2015 
# Modified 13 Apr 2016

##############################################################################
Sys.time()
##############################################################################

library(iCOBRA)
library(Hmisc)
library(DEXSeq)
library(DRIMSeq)
library(plyr)
library(limma)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/gosia/multinomial_project/simulations_sim5'
# simulation='hsapiens_node_nonull'
# count_method='kallistofiltered5'
# filter_method="filter0"
# CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R'
# method_out='drimseq_0_3_3'
# comparison_out='drimseq_0_3_3_comparison'


##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(simulation)
print(count_method)
print(filter_method)


##############################################################################

setwd(paste0(rwd, "/", simulation))



out_dir <- paste0(comparison_out, "/", filter_method, "/", count_method, "_")

dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE)

### colors

load(paste0(rwd, "/colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)


#######################################################
# Redo the dispersion plots with the same axes and larger font 
#######################################################


results_dir <- paste0(method_out, "/", count_method, "/", filter_method, "/")
files <- list.files(path = results_dir, pattern = "_d.Rdata" )
files


### For fitting smooth line, do not take into account the boudry grid points
disp_init <- common_disp <- as.numeric(read.table(paste0(results_dir, "common_dispersion.txt")))
disp_grid_length = 21 
disp_grid_range = c(-10, 10)
splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length)
splineDisp <- disp_init * 2^splinePts


for(i in grep("_grid_", files)){
  # i = 3
  
  load(paste0(results_dir, files[i]))
  
  text_size <- 26
  
  ggp <- plotDispersion(d) +
    coord_cartesian(ylim = log10(c(splineDisp[1], splineDisp[disp_grid_length]))) +
    theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = text_size), legend.title = element_text(size = text_size, face = "bold"))
  
  if(grepl("_grid_none_", files[i])){
    
    ggdf <- ggp$data
    not_boundry <- ggdf$dispersion < log10(splineDisp[disp_grid_length - 3]) & ggdf$dispersion > log10(splineDisp[1])
    ggdf <- ggdf[not_boundry, ]
    
    ggp <- ggp + 
      geom_smooth(data = ggdf, color = "black", span = 0.1)  
    
  }
  
  pdf(paste0(results_dir, gsub("_d.Rdata", "", files[i]), "_dispersion_vs_mean.pdf"))
  print(ggp)
  dev.off()
  
  ### Plot p-values
  ggp <- plotTest(d) +
    theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), plot.title = element_text(size = text_size))
  
  pdf(paste0(results_dir, gsub("_d.Rdata", "", files[i]), "_hist_pvalues.pdf"))
  print(ggp)
  dev.off()
  
  
}



#######################################################
# merge results for iCOBRA
#######################################################

####################### results from DEXSeq

if(count_method == "htseq")
  results_dir <- "4_results/dexseq_htseq_nomerge"
if(count_method == "kallisto")
  results_dir <- "4_results/dexseq_kallisto"
if(count_method == "htseqprefiltered15")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast15"
if(count_method == "htseqprefiltered5")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast5"
if(count_method == "kallistofiltered5")
  results_dir <- "4_results/dexseq_kallisto_txfilt_5"
if(count_method == "kallistoprefiltered5")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_kallisto_kallistoest_atleast5"

results_dir


### Make a histogram of DEXSeq exon p-values

# load(paste0(results_dir, ".Rdata"))
# 
# 
# ggp <- DRIMSeq:::dm_plotPvalues(pvalues = res$pvalue) +
#   geom_histogram(breaks = seq(0,1, by = 0.01), fill = "grey50")
# 
# pdf(paste0(out_dir, "dexseq_hist_pvalues_exons.pdf"))
# print(ggp)
# dev.off()



# pdf("test_dexseq_plot.pdf")
# 
# try(plotDEXSeq(res, geneID = "ENSG00000000003", FDR = 0.1, fitExpToVar = "condition", norCounts=FALSE, expression=TRUE, splicing = TRUE,  displayTranscripts=FALSE, names=FALSE, legend=TRUE, color=NULL, color.samples=NULL, las = 3), silent = TRUE)
# 
# dev.off()



results_padj <- list()


rt <- read.table(paste0(results_dir, ".txt"), header = TRUE, as.is = TRUE)
head(rt)

colnames(rt) <- c("gene_id", "dexseq")
head(rt)

results_padj[["dexseq"]] <- rt


####################### results from DRIMSeq

results_dir <- paste0(method_out, "/", count_method, "/", filter_method, "/")
files <- list.files(path = results_dir, pattern = "_results.txt" )
files

for(i in 1:length(files)){
  # i = 1
  method_name <- gsub(pattern = "_results.txt", replacement = "", x = files[i])
  
  rt <- read.table(paste0(results_dir, files[i]), header = TRUE, as.is = TRUE)
  head(rt)
  
  rt <- rt[,c("gene_id","adj_pvalue")]
  colnames(rt) <- c("gene_id", method_name)
  
  results_padj[[method_name]] <- rt 
  
}


results_padj <- Reduce(function(...) merge(..., by = "gene_id", all=TRUE, sort = FALSE), results_padj)
rownames(results_padj) <- results_padj$gene_id




results_padj <- results_padj[, colnames(results_padj) %in% colors_df$methods]

keep_methods <- colors_df$methods %in% colnames(results_padj)

clolors <- colors[keep_methods]
colors_df <- colors_df[keep_methods, , drop = FALSE]

results_padj <- results_padj[, colors_df$methods]



#######################################################
# load simulation info
#######################################################


# simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, as.is = TRUE, sep = "\t")

truth_file <- list.files("3_truth/", pattern = "truth")
truth_file

truth <- read.table(paste0("3_truth/", truth_file), header = TRUE, as.is = TRUE, sep = "\t")

truth <- truth[, c("gene", "ds_status", "de_status", "TPM", "nbr_isoforms", "diff_IsoPct", "nbrexonbins")]
rownames(truth) <- truth$gene


####################### stratification



### diff_IsoPct

truth$diff_IsoPct_cat <- cut2(truth$diff_IsoPct, cuts = c(1/4, 2/4))

n <- table(truth$diff_IsoPct_cat)
n
nds <- table(truth$diff_IsoPct_cat[truth$ds_status == 1])
nds

truth$diff_IsoPct_catn <- truth$diff_IsoPct_cat
levels(truth$diff_IsoPct_catn) <- paste0(levels(truth$diff_IsoPct_cat), " n = ", as.numeric(n), ", nds = ", as.numeric(nds))
levels(truth$diff_IsoPct_catn)


### nbr_isoforms

table(truth$nbr_isoforms)

cuts <- quantile(truth$nbr_isoforms[truth$ds_status == 1], probs = c(1/3, 2/3), na.rm = TRUE)
cuts

# cuts <- quantile(truth$nbr_isoforms[truth$ds_status == 1 & truth$nbr_isoforms > 3], probs = c(1/2), na.rm = TRUE)
# cuts <- c(4, cuts)
# cuts


truth$nbr_isoforms_cat <- cut2(truth$nbr_isoforms, cuts = cuts)

n <- table(truth$nbr_isoforms_cat)
n
nds <- table(truth$nbr_isoforms_cat[truth$ds_status == 1])
nds


truth$nbr_isoforms_catn <- truth$nbr_isoforms_cat
levels(truth$nbr_isoforms_catn) <- paste0(levels(truth$nbr_isoforms_cat), " n = ", as.numeric(n), ", nds = ", as.numeric(nds))
levels(truth$nbr_isoforms_catn)



### nbrexonbins

table(truth$nbrexonbins)

cuts <- quantile(truth$nbrexonbins[truth$ds_status == 1], probs = c(1/3, 2/3), na.rm = TRUE)
cuts


truth$nbrexonbins_cat <- cut2(truth$nbrexonbins, cuts = cuts)

n <- table(truth$nbrexonbins_cat)
n
nds <- table(truth$nbrexonbins_cat[truth$ds_status == 1])
nds

truth$nbrexonbins_catn <- truth$nbrexonbins_cat
levels(truth$nbrexonbins_catn) <- paste0(levels(truth$nbrexonbins_cat), " n = ", as.numeric(n), ", nds = ", as.numeric(nds))
levels(truth$nbrexonbins_catn)


##########################################################################
# Load counts
##########################################################################


if(count_method == "htseq")
  count_dir <- "2_counts/dexseq_nomerge/dexseq"
if(count_method == "kallisto")
  count_dir <- "2_counts/kallisto/kallisto"
if(count_method == "htseqprefiltered15")
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast15/dexseq"
if(count_method == "htseqprefiltered5")
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast5/dexseq"
if(count_method == "kallistofiltered5")
  count_dir <- "2_counts/kallisto_txfilt_5/kallisto"
if(count_method == "kallistoprefiltered5")
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/kallisto_kallistoest_atleast5/kallisto"

count_dir

### load counts
counts_list <- lapply(1:6, function(i){
  # i = 1
  cts <- read.table(paste0(count_dir, i, ".txt"), header = FALSE, as.is = TRUE)
  colnames(cts) <- c("group_id", paste0("sample_", i))  
  return(cts)
})

counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)

counts <- counts[!grepl(pattern = "_", counts$group_id),]
group_split <- strsplit2(counts[,1], ":")
counts <- counts[, -1]



#######################################################
# plots with iCOBRA
#######################################################

cobradata <- COBRAData(padj = results_padj, truth = truth)

save(cobradata, file = paste0(out_dir, "cobradata.Rdata"))




### Venn diagrams overall

cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", aspects = c("overlap"), thr_venn = 0.05)

basemethods(cobraperf)


for(i in grep("drimseq", basemethods(cobraperf))){
  # i = 2
  
  cobraplot <- prepare_data_for_plot(cobraperf, keepmethods = c("dexseq", basemethods(cobraperf)[i]), colorscheme = c(colors[c("dexseq", basemethods(cobraperf)[i])], "grey"), incltruth = TRUE)
  
  pdf(paste0(out_dir, "venn_", basemethods(cobraperf)[i], ".pdf"))
  plot_overlap(cobraplot, cex = c(1.2, 1, 0.7))
  dev.off()
  
  pdf(paste0(out_dir, "upset_", basemethods(cobraperf)[i], ".pdf"), width = 10)
  plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 14, sets = c("truth", "dexseq", basemethods(cobraperf)[i]), sets.bar.color = c("grey", colors[c("dexseq", basemethods(cobraperf)[i])]))
  dev.off()
  
  
}


### Plot distributions of differenet gene characteristics for unique genes
bmethds <- basemethods(cobraperf)

text_size=18 
legend_size=16


### Mean expression

# mean_expression <- truth$TPM
# names(mean_expression) <- truth$gene

mean_expression <- by(counts, factor(group_split[, 1]), function(x){
  mean(colSums(x, na.rm = TRUE), na.rm = TRUE)
}, simplify = FALSE)

mean_expression <- unlist(mean_expression)


### Number of expressed transcripts per gene 

# nbr_isoforms <- truth$nbr_isoforms
# names(nbr_isoforms) <- truth$gene


nbr_features<- by(counts, factor(group_split[, 1]), function(x){
  x <- as.matrix(x)
  sum(rowSums(x > 10, na.rm = TRUE) > 2, na.rm = TRUE)
}, simplify = FALSE)

nbr_features <- unlist(nbr_features)



for(i in grep("drimseq", bmethds)){
  # i = 2
  
  all_genes <- rownames(results_padj)
  truth_genes <- truth[truth$ds_status == 1, "gene"]
  
  m1 <- "dexseq"
  m2 <- bmethds[i]
  
  genes_sign_m1 <- rownames(results_padj[results_padj[, m1] < 0.05 & !is.na(results_padj[, m1]), , drop = FALSE])
  genes_sign_m2 <- rownames(results_padj[results_padj[, m2] < 0.05 & !is.na(results_padj[, m2]), , drop = FALSE])
  
  genes_sign_overlap <- intersect(genes_sign_m1, genes_sign_m2)
  genes_sign_m1_unique <- setdiff(genes_sign_m1, genes_sign_overlap)
  genes_sign_m2_unique <- setdiff(genes_sign_m2, genes_sign_overlap)
  
  
  ggdf <- data.frame(mean_expression = c(mean_expression[all_genes], mean_expression[truth_genes], mean_expression[genes_sign_m1], mean_expression[genes_sign_m2]), 
    group = c(rep("all_genes", length(all_genes)), rep("truth_genes", length(truth_genes)), rep(m1, length(genes_sign_m1)), rep(m2, length(genes_sign_m2))), 
    stringsAsFactors = FALSE)
  
  ggdf <- ggdf[complete.cases(ggdf), ]
  
  ggdf$group <- factor(ggdf$group, levels = c("all_genes", "truth_genes", m1, m2), , labels = paste0(c("all_genes", "truth_genes", m1, m2), " (", c(length(all_genes), length(truth_genes), length(genes_sign_m1), length(genes_sign_m2)), ")"))
  
  ggdf <- ggdf[ggdf$mean_expression > 0, ,drop = FALSE]
  
  
  ggp <- ggplot(ggdf, aes(x = log10(mean_expression), color = group, group = group)) +
    geom_density(size = 2) +
    theme_bw() +
    xlab("Log10 of gene mean expression") +
    theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = legend_size), legend.title = element_blank(), legend.position = "bottom") +
    scale_color_manual(values = as.character(c("grey50", "grey", colors[c(m1, m2)]))) +
    guides(color = guide_legend(nrow = 2))
  
  pdf(paste0(out_dir, "characteristics_mean_expr_", bmethds[i], ".pdf"))
  print(ggp)
  dev.off()
  

  ggdf <- data.frame(nbr_features = nbr_features[c(all_genes, truth_genes, genes_sign_m1, genes_sign_m2)], 
    group = rep(c("all_genes", "truth_genes", m1, m2), c(length(all_genes), length(truth_genes), length(genes_sign_m1), length(genes_sign_m2))),
    stringsAsFactors = FALSE)
  
  ggdf <- ggdf[complete.cases(ggdf), ]
  
  ggdf$group <- factor(ggdf$group, levels = c("all_genes", "truth_genes", m1, m2), , labels = paste0(c("all_genes", "truth_genes", m1, m2), " (", c(length(all_genes), length(truth_genes), length(genes_sign_m1), length(genes_sign_m2)), ")"))
  
  
  ggp <- ggplot(ggdf, aes(x = nbr_features, color = group, group = group)) +
    geom_density(size = 2) +
    theme_bw() +
    xlab("Number of expressed features") +
    theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = legend_size), legend.title = element_blank(), legend.position = "bottom") +
    scale_color_manual(values = as.character(c("grey50", "grey", colors[c(m1, m2)])))+
    guides(color = guide_legend(nrow = 2))
  
  pdf(paste0(out_dir, "characteristics_nbr_features_", bmethds[i], ".pdf"))
  print(ggp)
  dev.off()
  
  
}




### FDR TPR overall

cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", aspects = c("fdrtpr", "fdrtprcurve"), onlyshared = FALSE, maxsplit = Inf)

cobraplot <- prepare_data_for_plot(cobraperf, keepmethods = factor(colors_df$methods, levels = colors_df$methods), colorscheme = colors[basemethods(cobraperf)])

cobraplot <- iCOBRA:::reorder_levels(cobraplot, colors_df$methods)





## Plot only points
# xaxisrange <- range(cobraplot@fdrtpr$FDR, na.rm = TRUE)
# yaxisrange <- range(cobraplot@fdrtpr$TPR, na.rm = TRUE)
xaxisrange <- c(0, 0.7)
yaxisrange <- c(0.3, 1)

ggp <- plot_fdrtprcurve(cobraplot, plottype = c("points"), pointsize = 3, xaxisrange = c(0, xaxisrange[2]), yaxisrange = yaxisrange)
ggp <- ggp + 
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 10), strip.text = element_text(size = 11)) + 
  guides(colour = guide_legend(nrow = 2)) 
  

pdf(paste0(out_dir, "fdrtpr.pdf"), 7, 8)
print(ggp)
dev.off()



## Plot the points and the curve
# xaxisrange <- range(cobraplot@fdrtprcurve$FDR, na.rm = TRUE)
# yaxisrange <- range(cobraplot@fdrtprcurve$TPR[cobraplot@fdrtprcurve$TPR > 0], na.rm = TRUE)
xaxisrange <- c(0, 0.8)
yaxisrange <- c(0.2, 1)

ggp <- plot_fdrtprcurve(cobraplot, plottype = c("curve", "points"), pointsize = 3, xaxisrange = c(0, xaxisrange[2]), yaxisrange = yaxisrange)
ggp <- ggp + 
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 10), strip.text = element_text(size = 11)) + 
  guides(colour = guide_legend(nrow = 2))

pdf(paste0(out_dir, "fdrtprcurve.pdf"), 7, 8)
print(ggp)
dev.off()



### FDR TPR stratified

splv_list <- c("diff_IsoPct_catn", "nbr_isoforms_catn", "nbrexonbins_catn")

for(i in 1:length(splv_list)){
  # i = 3
  
  splv <- splv_list[i]
  
  cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", splv = splv, c("fdrtpr", "fdrtprcurve"), onlyshared = FALSE, maxsplit = Inf)
  
  cobraplot <- prepare_data_for_plot(cobraperf, incloverall = TRUE, colorscheme = colors[basemethods(cobraperf)])
  
  cobraplot <- iCOBRA:::reorder_levels(cobraplot, colors_df$methods)
  
  ### Plot only points
  # levels(cobraplot@fdrtpr$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtpr$splitval))
  levels(cobraplot@fdrtpr$splitval) <- gsub(paste0("_catn"), "", levels(cobraplot@fdrtpr$splitval))
  
  xaxisrange <- range(cobraplot@fdrtpr$FDR, na.rm = TRUE)
  yaxisrange <- range(cobraplot@fdrtpr$TPR, na.rm = TRUE)
  
  ggp <- plot_fdrtprcurve(cobraplot, plottype = c("points"), pointsize = 3, stripsize = 8, xaxisrange = c(0, xaxisrange[2]), yaxisrange = yaxisrange)
  ggp <- ggp + 
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), strip.text = element_text(size = 12)) + 
    guides(colour = guide_legend(nrow = 1)) + 
    facet_wrap(~splitval, nrow = 1) 
  
  pdf(paste0(out_dir, "fdrtpr_", splv ,"_auto.pdf"), 14, 5)
  print(ggp)
  dev.off()
  
  
  xaxisrange <- c(0, 0.6)
  yaxisrange <- c(0.3, 0.9)
  
  ggp <- plot_fdrtprcurve(cobraplot, plottype = c("points"), pointsize = 3, stripsize = 8, xaxisrange = c(0, xaxisrange[2]), yaxisrange = yaxisrange)
  ggp <- ggp + 
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), strip.text = element_text(size = 12)) + 
    guides(colour = guide_legend(nrow = 1)) + 
    facet_wrap(~splitval, nrow = 1) 
  
  pdf(paste0(out_dir, "fdrtpr_", splv ,".pdf"), 14, 5)
  print(ggp)
  dev.off()
  
  
  
  
  ### Plot the points and the curve
  # levels(cobraplot@fdrnbrcurve$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrnbrcurve$splitval))
  levels(cobraplot@fdrnbrcurve$splitval) <- gsub(paste0("_catn"), "", levels(cobraplot@fdrnbrcurve$splitval))
  # levels(cobraplot@fdrtprcurve$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtprcurve$splitval))
  levels(cobraplot@fdrtprcurve$splitval) <- gsub(paste0("_catn"), "", levels(cobraplot@fdrtprcurve$splitval))
  
  xaxisrange <- range(cobraplot@fdrtprcurve$FDR, na.rm = TRUE)
  yaxisrange <- range(cobraplot@fdrtprcurve$TPR[cobraplot@fdrtprcurve$TPR > 0], na.rm = TRUE)

  ggp <- plot_fdrtprcurve(cobraplot, plottype = c("curve", "points"), pointsize = 3, stripsize = 8, xaxisrange = c(0, xaxisrange[2]), yaxisrange = yaxisrange)
  ggp <- ggp + 
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), strip.text = element_text(size = 12)) + 
    guides(colour = guide_legend(nrow = 1)) + 
    facet_wrap(~splitval, nrow = 1)
  
  pdf(paste0(out_dir, "fdrtprcurve_", splv ,"_auto.pdf"), 14, 5)
  print(ggp)
  dev.off()
  
  
  xaxisrange <- c(0, 0.6)
  yaxisrange <- c(0.3, 0.9)
  
  ggp <- plot_fdrtprcurve(cobraplot, plottype = c("curve", "points"), pointsize = 3, stripsize = 8, xaxisrange = c(0, xaxisrange[2]), yaxisrange = yaxisrange)
  ggp <- ggp + 
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), strip.text = element_text(size = 12)) + 
    guides(colour = guide_legend(nrow = 1)) + 
    facet_wrap(~splitval, nrow = 1) 
  
  pdf(paste0(out_dir, "fdrtprcurve_", splv ,".pdf"), 14, 5)
  print(ggp)
  dev.off()
  
  
  
}


#######################################################
# CAT plot
#######################################################


results_padj <- list()
metadata <- list()


### results from DEXSeq
if(count_method == "htseq")
  results_dir <- "4_results/dexseq_htseq_nomerge"
if(count_method == "kallisto")
  results_dir <- "4_results/dexseq_kallisto"
if(count_method == "htseqprefiltered15")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast15"
if(count_method == "htseqprefiltered5")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast5"
if(count_method == "kallistofiltered5")
  results_dir <- "4_results/dexseq_kallisto_txfilt_5"
if(count_method == "kallistoprefiltered5")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_kallisto_kallistoest_atleast5"

results_dir


rt <- read.table(paste0(results_dir, ".txt"), header = TRUE, as.is = TRUE)
head(rt)

colnames(rt) <- c("gene_id", "dexseq")
head(rt)


method_name <- "dexseq"
results_padj[[method_name]] <- rt
metadata[[method_name]] <- data.frame(method_name = method_name, stringsAsFactors = FALSE)


### results from DRIMSeq

results_dir <- paste0(method_out, "/", count_method, "/", filter_method, "/")
files <- list.files(path = results_dir, pattern = "_results.txt" )
files

for(i in 1:length(files)){
  # i = 1
  method_name <- gsub(pattern = "_results.txt", replacement = "", x = files[i])
  
  rt <- read.table(paste0(results_dir, files[i]), header = TRUE, as.is = TRUE)
  head(rt)
  
  rt <- rt[,c("gene_id","adj_pvalue")]
  colnames(rt) <- c("gene_id", method_name)
  
  results_padj[[method_name]] <- rt 
  metadata[[method_name]] <- data.frame(method_name = method_name, stringsAsFactors = FALSE)
  
}


metadata <- rbind.fill(metadata)

metadata$method_name <- factor(metadata$method_name, levels = colors_df$methods)
metadata$method_name <- factor(metadata$method_name)
metadata


names(results_padj)

### Create results list using adjusted p-values as p-values
results_padj <- lapply(results_padj, function(x){
  
  colnames(x) <- c("gene_id", "pvalue")
  x$adj_pvalue <- x$pvalue
  
  return(x)
  
})




source(CAT_function_path)

### Index that indicates pairs of methods to be compared
metadata$results_order <- 1:nrow(metadata)

reference_method <- "dexseq"

indx <- lapply(1, function(i){
  # i = 1
  
  metadata_tmp <- metadata
  
  indx1 <- metadata_tmp[metadata_tmp$method_name == reference_method, "results_order"]
  indx2 <- metadata_tmp[!metadata_tmp$method_name == reference_method, "results_order"]
  
  data.frame(indx1 = rep(indx1, length(indx2)), indx2 = indx2)
  
})

indx <- rbind.fill(indx)
indx


### Calculate CAT data per pair of methods
data_CAT <- lapply(1:nrow(indx), function(j){
  # j = 4
  
  calculateCAT(results1 = results_padj[[indx$indx1[j]]], results2 = results_padj[[indx$indx2[j]]], by = 10)
  
})


### Metadata for overlaps in data_Overlaps list
metadata_ov <- metadata[indx$indx2, ]
metadata_ov$method_name <- factor(metadata_ov$method_name)


ggp <- plotCAT(data_CAT, metadata = metadata_ov, plot_var = "method_name", facet_var = NULL, plot_colors = colors[levels(metadata_ov$method_name)], plotx = TRUE, reference_color = colors[reference_method])

ggp <- ggp + 
  coord_cartesian(xlim = c(0, 2000), ylim = c(0, 1)) 


pdf(paste0(out_dir, "cat.pdf"), 7, 7)
print(ggp)
dev.off()















































