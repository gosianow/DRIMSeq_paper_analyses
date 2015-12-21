######################################################
## ----- sim5_drimseq_comparison_filtering
## <<sim5_drimseq_comparison_filtering.R>>

# BioC 3.2
# Created 18 Dec 2015 


##############################################################################

library(iCOBRA)
library(Hmisc)
library(plyr)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/gosia/multinomial_project/simulations_sim5'
# 
# simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull','hsapiens_withde_nonull')
# count_method_list=c('kallisto','htseq')
# filter_method_list=c('filter0','filter1','filter2','filter3')
# name=''
# legend_nrow=1
# pdf_width=12
# pdf_height=7

##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

print(rwd)
print(simulation_list)
print(count_method_list)
print(filter_method_list)
print(name)
print(legend_nrow)
print(pdf_width)
print(pdf_height)


##############################################################################

setwd(rwd)
method_out <- "drimseq_0_3_3"


### colors

load(paste0(rwd, "/colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)

### Plot

out_dir_plots <- paste0("drimseq_0_3_3_comparison/", "different_filtering", "/")
dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)

##############################################################################

### Load results

results_padj_list <- list()
truth_list <- list()
ix <- 1

for(i in 1:length(simulation_list)){
  
  for(j in 1:length(count_method_list)){
    # i = 1; j = 1;
    
    results_padj_tmp <- list()
    
    for(k in 1:length(filter_method_list)){
      # i = 1; j = 1; k = 3
      
      simulation <- simulation_list[i]
      count_method <- count_method_list[j]
      filter_method <- filter_method_list[k]
      
      comparison_out <- paste0(simulation, "/drimseq_0_3_3_comparison")
      out_dir <- paste0(comparison_out, "/", filter_method, "/", count_method, "_")
      
      if(!file.exists(paste0(out_dir, "cobradata.Rdata")))
        next
      
      load(paste0(out_dir, "cobradata.Rdata"))
      
      results_padj <- cobradata@padj
      colnames(results_padj) <- paste0(colnames(results_padj), ".", filter_method)
      results_padj$rns <- paste0(simulation, "_", count_method, ".", rownames(cobradata@padj))
      results_padj_tmp[[paste0(ix)]] <- results_padj
      
      if(k == 1){
        truth <-  cobradata@truth
        rownames(truth) <- paste0(simulation, "_", count_method, ".", rownames(cobradata@truth))
        truth$simulation <- simulation
        truth$count_method <- count_method
        truth_list[[paste0(ix)]] <- truth
      }
      
      ix <- ix + 1
    }
    
    results_padj_tmp <- Reduce(function(...) merge(..., by = "rns", all = TRUE, sort = FALSE), results_padj_tmp)
    rownames(results_padj_tmp) <- results_padj_tmp$rns
    
    results_padj_list[[ix]] <- results_padj_tmp[, -1]
    
  }
}




results_padj <- rbind.fill(results_padj_list)
rownames(results_padj) <- unlist(lapply(results_padj_list, rownames))

truth <- rbind.fill(truth_list)
rownames(truth) <- unlist(lapply(truth_list, rownames))


table(truth$simulation)
table(truth$count_method)


truth$simulation <- factor(truth$simulation, levels = simulation_list)
truth$count_method <- factor(truth$count_method, levels = count_method_list)

all(rownames(results_padj) %in% rownames(truth))


results_padj <- results_padj[, -grep(pattern = "dexseq", colnames(results_padj))[-1]]

colnames(results_padj)[grep(pattern = "dexseq", colnames(results_padj))] <- "dexseq"



##############################################################################
### Order for split - simulations per row
##############################################################################


truth$split <- interaction(truth$simulation, truth$count_method, lex.order = FALSE)

levels(truth$split)

facet_nrow <- length(count_method_list)




for(drimseq_version in c("drimseq_genewise_grid_none", "drimseq_genewise_grid_common")){
  # drimseq_version <- "drimseq_genewise_grid_none"

  
  cobradata <- COBRAData(padj = results_padj[, grep(paste0("dexseq|", drimseq_version), colnames(results_padj))], truth = truth)
  
  
  ### FDR TPR stratified
  
  splv <- "split"
  
  cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", splv = splv, aspects = c("fdrtpr", "fdrtprcurve"), onlyshared = FALSE, maxsplit = Inf)
  
  cobraplot <- prepare_data_for_plot(cobraperf, incloverall = FALSE)
  
  
  levels(cobraplot@fdrtpr$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtpr$splitval))
  
  
  ggp <- plot_fdrtprcurve(cobraplot, plottype = c("points"), pointsize = 3, stripsize = 9, xaxisrange = c(0, 0.7), yaxisrange = c(0.4, 1))
  ggp <- ggp + 
    theme(legend.position = "bottom", strip.text = element_text(size = 11)) + 
    guides(colour = guide_legend(nrow = legend_nrow)) + 
    facet_wrap(~splitval, nrow = facet_nrow)
  
  
  pdf(paste0(out_dir_plots, "fdrtpr_", drimseq_version , name, ".pdf"), width = pdf_width, height = pdf_height)
  print(ggp)
  dev.off()
  
  
  
  levels(cobraplot@fdrnbrcurve$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrnbrcurve$splitval))
  levels(cobraplot@fdrtprcurve$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtprcurve$splitval))
  
  
  ggp <- plot_fdrtprcurve(cobraplot, pointsize = 3, stripsize = 9, xaxisrange = c(0, 1), yaxisrange = c(0, 1))
  ggp <- ggp + 
    theme(legend.position = "bottom", strip.text = element_text(size = 11)) + 
    guides(colour = guide_legend(nrow = legend_nrow)) + 
    facet_wrap(~splitval, nrow = facet_nrow)
  
  
  pdf(paste0(out_dir_plots, "fdrtprcurve_", drimseq_version, name, ".pdf"), width = pdf_width, height = pdf_height)
  print(ggp)
  dev.off()
  
}










