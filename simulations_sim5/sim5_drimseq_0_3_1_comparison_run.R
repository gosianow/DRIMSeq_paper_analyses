######################################################
## ----- sim5_drimseq_0_3_1_comparison_run
## <<sim5_drimseq_0_3_1_comparison_run.R>>

# BioC 3.1
# Created 16 Nov 2015 

##############################################################################

library(iCOBRA)
library(Hmisc)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_sim5'
simulation='hsapiens_node_nonull'
count_method=c('htseq','kallisto','htseq_prefiltered15')[1]


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



##############################################################################


setwd(paste0(rwd, "/", simulation))
method_out <- "drimseq_0_3_1"

comparison_out <- "drimseq_0_3_1_comparison"

out_dir <- paste0(comparison_out, "/", count_method, "_")

dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE)

### colors

load(paste0(rwd, "/colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)

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




results_padj <- list()


rt <- read.table(paste0(results_dir, ".txt"), header = TRUE, as.is = TRUE)
head(rt)

colnames(rt) <- c("gene_id", "dexseq")
head(rt)

results_padj[["dexseq"]] <- rt


####################### results from DRIMSeq

results_dir <- paste0(method_out, "/", count_method, "/")
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


truth$nbr_isoforms_cat <- cut2(truth$nbr_isoforms, cuts = cuts)

n <- table(truth$nbr_isoforms_cat)
n
nds <- table(truth$nbr_isoforms_cat[truth$ds_status == 1])
nds


truth$nbr_isoforms_catn <- truth$nbr_isoforms_cat
levels(truth$nbr_isoforms_catn) <- paste0(levels(truth$nbr_isoforms_cat), " n = ", as.numeric(n), ", nds = ", as.numeric(nds))
levels(truth$nbr_isoforms_catn)



### nbr_isoforms

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





#######################################################
# plot with iCOBRA
#######################################################

cobradata <- COBRAData(padj = results_padj, truth = truth)

save(cobradata, file = paste0(out_dir, "cobradata.Rdata"))




### Venn diagrams overall

cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", aspects = c("overlap"))

basemethods(cobraperf)


for(i in 2:(length(basemethods(cobraperf)) -1)){
  # i = 2
  
  cobraplot <- prepare_data_for_plot(cobraperf, keepmethods = c("dexseq", basemethods(cobraperf)[i]), colorscheme = c(colors[c("dexseq", basemethods(cobraperf)[i])], "grey"), incltruth = TRUE)
  
  pdf(paste0(out_dir, "venn_", basemethods(cobraperf)[i], ".pdf"))
  plot_overlap(cobraplot, cex = c(1.2, 1, 0.7))
  dev.off()
  
}



### FDR TPR overall

cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", aspects = "fdrtpr")


cobraplot <- prepare_data_for_plot(cobraperf, keepmethods = factor(colors_df$methods, levels = colors_df$methods), colorscheme = colors[basemethods(cobraperf)])


cobraplot@fdrtpr$method <- factor(cobraplot@fdrtpr$method, levels = colors_df$methods)


ggp <- plot_fdrtprcurve(cobraplot, plottype = c("points"), pointsize = 3, xaxisrange = c(0, 0.6), yaxisrange = c(0.4, 1))
ggp <- ggp + 
  theme(legend.position = "right")

pdf(paste0(out_dir, "fdrtpr.pdf"), 7, 5)
print(ggp)
dev.off()




### FDR TPR stratified

splv_list <- c("diff_IsoPct_catn", "nbr_isoforms_catn", "nbrexonbins_catn")

for(i in 1:length(splv_list)){
  # i = 3
  
  splv <- splv_list[i]
  
  cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", splv = splv, aspects = "fdrtpr", onlyshared = FALSE, maxsplit = Inf)
  
  cobraplot <- prepare_data_for_plot(cobraperf, incloverall = FALSE, colorscheme = colors[basemethods(cobraperf)])
  
  cobraplot@fdrtpr$method <- factor(cobraplot@fdrtpr$method, levels = colors_df$methods)
  
  
  # levels(cobraplot@fdrtpr$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtpr$splitval))
  levels(cobraplot@fdrtpr$splitval) <- gsub(paste0("_catn"), "", levels(cobraplot@fdrtpr$splitval))
  
  ggp <- plot_fdrtprcurve(cobraplot, plottype = c("points"), pointsize = 3, stripsize = 8, xaxisrange = c(0, 0.7), yaxisrange = c(0.3, 1))
  ggp <- ggp + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(nrow = 2)) + 
    facet_wrap(~splitval, nrow = 1)
  
  pdf(paste0(out_dir, "fdrtpr_", splv ,".pdf"), 10, 5)
  print(ggp)
  dev.off()
  
}



# save(cobradata, colors, file = "/home/gosia/case_for_Charlotte.Rdata")
# 
# levels(truth$diff_IsoPct_catn)
# 
# levels(cobradata@truth$diff_IsoPct_catn)





