######################################################
## ----- sim5_drimseq_comparison_filtering
## <<sim5_drimseq_comparison_filtering.R>>

# BioC 3.2
# Created 18 Dec 2015 
# Modified 13 Apr 2016

##############################################################################

library(ggplot2)
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
# prefilter_method_list=c('kallistoprefiltered5','htseqprefiltered5')
# name='_prefilt'
# legend_nrow=2
# pdf_width=9
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
      
      rpadj <- cobradata@padj
      colnames(rpadj) <- paste0(colnames(rpadj), ".", filter_method)
      rpadj$rns <- paste0(simulation, "_", count_method, ".", rownames(cobradata@padj))
      results_padj_tmp[[paste0(ix)]] <- rpadj
      
      if(k == 1){
        tr <-  cobradata@truth
        rownames(tr) <- paste0(simulation, "_", count_method, ".", rownames(cobradata@truth))
        tr$simulation <- simulation
        tr$count_method <- count_method
        truth_list[[paste0(ix)]] <- tr
      }
      
      ix <- ix + 1
    }
    
    results_padj_tmp <- Reduce(function(...) merge(..., by = "rns", all = TRUE, sort = FALSE), results_padj_tmp)
    rownames(results_padj_tmp) <- results_padj_tmp$rns
    
    results_padj_list[[ix]] <- results_padj_tmp[, -grep("rns", colnames(results_padj_tmp))]
    
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


results_padj <- results_padj[, -grep(pattern = "dexseq.filter[1-9]", colnames(results_padj))]




### add results of prefiltering with filter0
if(!is.null(prefilter_method_list)){
  
  count_method_list <- prefilter_method_list
  filter_method_list <- "filter0"
  simulation_list <- simulation_list[!grepl("hsapiens_withde", simulation_list)]
    
  results_padj_list <- list()
  truth_list <- list()
  ix <- 1
  
  for(i in 1:length(simulation_list)){
    
    for(j in 1:length(count_method_list)){
      # i = 1; j = 1;
      
      results_padj_tmp <- list()
      
      for(k in 1:length(filter_method_list)){
        # i = 2; j = 2; k = 1
        
        simulation <- simulation_list[i]
        count_method <- count_method_list[j]
        filter_method <- filter_method_list[k]
        
        comparison_out <- paste0(simulation, "/drimseq_0_3_3_comparison")
        out_dir <- paste0(comparison_out, "/", filter_method, "/", count_method, "_")
        
        if(!file.exists(paste0(out_dir, "cobradata.Rdata")))
          next
        
        load(paste0(out_dir, "cobradata.Rdata"))
        
        rpadj <- cobradata@padj
        colnames(rpadj) <- paste0(colnames(rpadj), ".", filter_method)
        rpadj$rns <- paste0(simulation, "_", count_method, ".", rownames(cobradata@padj))
        results_padj_tmp[[paste0(ix)]] <- rpadj
        
        print(out_dir)
        print(head(rpadj))
        
        if(k == 1){
          tr <-  cobradata@truth
          rownames(tr) <- paste0(simulation, "_", count_method, ".", rownames(cobradata@truth))
          tr$simulation <- simulation
          tr$count_method <- count_method
          truth_list[[paste0(ix)]] <- tr
        }
        
        ix <- ix + 1
      }
      
      results_padj_tmp <- Reduce(function(...) merge(..., by = "rns", all = TRUE, sort = FALSE), results_padj_tmp)
      rownames(results_padj_tmp) <- results_padj_tmp$rns
      
      results_padj_list[[ix]] <- results_padj_tmp[, -grep("rns", colnames(results_padj_tmp))]
      
      print(head(results_padj_list[[ix]]))
      
    }
  }
  
  
  pre_results_padj <- rbind.fill(results_padj_list)
  rownames(pre_results_padj) <- unlist(lapply(results_padj_list, rownames))
  
  rownames(pre_results_padj) <- gsub("prefiltered5", "", rownames(pre_results_padj))
  
  colnames(pre_results_padj) <- gsub("filter0", "prefilter5", colnames(pre_results_padj))
  
  results_padj <- merge(results_padj, pre_results_padj, by = 0, all = TRUE, sort = FALSE)
  
  rownames(results_padj) <- results_padj$Row.names
  
  results_padj <- results_padj[, -grep("Row.names", colnames(results_padj))]
  
}



all(rownames(results_padj) %in% rownames(truth))





##############################################################################
### Order for split - simulations per row
##############################################################################

truth <- truth[rownames(results_padj), ]


truth$split <- interaction(truth$simulation, truth$count_method, lex.order = FALSE)

levels(truth$split)

head(results_padj[truth$split == "hsapiens_node_nonull.htseq", "dexseq.prefilter5"])



facet_nrow <- length(count_method_list)




for(drimseq_version in c("drimseq_genewise_grid_none", "drimseq_genewise_grid_common", "drimseq_genewise_grid_trended")){
  # drimseq_version <- "drimseq_genewise_grid_none"
  
  
  cobradata <- COBRAData(padj = results_padj[, grep(paste0("dexseq|", drimseq_version), colnames(results_padj))], truth = truth)
  
  
  ### FDR TPR stratified
  
  splv <- "split"
  
  cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", splv = splv, aspects = c("fdrtpr", "fdrtprcurve"), onlyshared = FALSE, maxsplit = Inf)
  
  cobraplot <- prepare_data_for_plot(cobraperf, incloverall = FALSE)
  
  
  
  ### Plot points only
  
  levels(cobraplot@fdrtpr$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtpr$splitval))
  
  xaxisrange = c(0, 0.7)
  yaxisrange = c(0.4, 0.8)
  
  ggp <- plot_fdrtprcurve(cobraplot, plottype = c("points"), pointsize = 3)
  ggp <- ggp + 
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16, color = "darkgrey"), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 10), strip.text = element_text(size = 11), strip.background = element_rect(colour = "black", fill="white")) + 
    guides(colour = guide_legend(nrow = legend_nrow)) + 
    facet_wrap(~splitval, nrow = facet_nrow) +
    coord_cartesian(xlim = xaxisrange, ylim = yaxisrange)
  
  
  pdf(paste0(out_dir_plots, "fdrtpr_", drimseq_version , name, ".pdf"), width = pdf_width, height = pdf_height)
  print(ggp)
  dev.off()
  
  
  ### Plot points and a curve
  
  levels(cobraplot@fdrnbrcurve$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrnbrcurve$splitval))
  levels(cobraplot@fdrtprcurve$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtprcurve$splitval))
  
  xaxisrange = c(0, 0.7)
  yaxisrange = c(0.4, 0.8)
  
  ggp <- plot_fdrtprcurve(cobraplot, pointsize = 3)
  ggp <- ggp + 
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16, color = "darkgrey"), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 10), strip.text = element_text(size = 11), strip.background = element_rect(colour = "black", fill="white")) + 
    guides(colour = guide_legend(nrow = legend_nrow)) + 
    facet_wrap(~splitval, nrow = facet_nrow) +
    coord_cartesian(xlim = xaxisrange, ylim = yaxisrange)
  
  
  pdf(paste0(out_dir_plots, "fdrtprcurve_", drimseq_version, name, ".pdf"), width = pdf_width, height = pdf_height)
  print(ggp)
  dev.off()
  
}










