######################################################
## ----- kim_drimseq_0_3_1_comparison_run
## <<kim_drimseq_0_3_1_comparison_run.R>>

# BioC 3.1
# Created 9 Nov 2015 

##############################################################################

library(DRIMSeq)
library(iCOBRA)

source("/home/gosia/R/drimseq_paper/dm_comparison_functions/dm_plotVenn.R")

##############################################################################
# Read in the arguments
##############################################################################

# rwd='/home/Shared/data/seq/kim_adenocarcinoma/'
# count_method=c('htseq', 'kallisto')[1]
# model=c('model_full', 'model_null_normal1', 'model_null_tumor1')[1]


## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(model)
print(count_method)


##############################################################################

setwd(rwd)
method_out <- "drimseq_0_3_1"
comparison_out <- "drimseq_0_3_1_comparison"


out_dir <- paste0(comparison_out, "/",  model, "/", count_method, "/")
dir.create(out_dir, recursive = TRUE)



#######################################################
# merge results for iCOBRA
#######################################################


results_padj <- list()


####################### results from DEXSeq

rt <- read.table(paste0("4_results/dexseq_1_10_8/", model,"/", count_method, "/dexseq_gene_results.txt"), header = TRUE, as.is = TRUE)
head(rt)

colnames(rt) <- c("gene_id", "dexseq")
head(rt)

results_padj[["dexseq"]] <- rt


####################### results from DRIMSeq

res_path <- paste0(method_out, "/",  model, "/", count_method, "/")
files <- list.files(path = res_path, pattern = "_results.txt" )

for(i in 1:length(files)){
  # i = 1
  method_name <- gsub(pattern = "_results.txt", replacement = "", x = files[i])
  
  rt <- read.table(paste0(res_path, files[i]), header = TRUE, as.is = TRUE)
  head(rt)
  
  rt <- rt[,c("gene_id","adj_pvalue")]
  colnames(rt) <- c("gene_id", method_name)

  results_padj[[method_name]] <- rt 
 
}


results_padj <- Reduce(function(...) merge(..., by = "gene_id", all=TRUE, sort = FALSE), results_padj)
rownames(results_padj) <- results_padj$gene_id


####################### use iCOBRA


cobradata <- COBRAData(padj = results_padj[, -1])

cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 0.05)

basemethods(cobraperf)


for(i in 2:length(basemethods(cobraperf))){
  # i = 2
  
  cobraplot <- prepare_data_for_plot(cobraperf, keepmethods = c("dexseq", basemethods(cobraperf)[i]), colorscheme = c("dodgerblue3", "darkorange1"), incltruth = FALSE)
  
  pdf(paste0(out_dir, "/venn_", basemethods(cobraperf)[i], ".pdf"))
  plot_overlap(cobraplot, cex=c(1.2,1,0.7))
  dev.off()
  
}





















































