######################################################
## ----- sim5_drimseq_comparison_combined
## <<sim5_drimseq_comparison_combined.R>>

# BioC 3.2
# Created 16 Nov 2015 
# Modified 13 Apr 2016

##############################################################################

library(iCOBRA)
library(Hmisc)
library(DEXSeq)
library(DRIMSeq)
library(plyr)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/gosia/multinomial_project/simulations_sim5'
# simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull')
# count_method_list=c('kallisto','kallistofiltered5','htseq','htseqprefiltered5')
# filter_method="filter0"
# CAT_function_path='/home/gosia/R/drimseq_paper/help_functions/dm_plotCAT.R'
# name=''
# legend_nrow=1
# pdf_width=14
# pdf_height=8
# method_out='drimseq_0_3_3'
# comparison_out='drimseq_0_3_3_comparison'

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
print(filter_method)
print(legend_nrow)
print(pdf_width)
print(pdf_height)


##############################################################################

setwd(rwd)



### colors

load(paste0(rwd, "/colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)

### Plot

out_dir_plots <- paste0(comparison_out, "/", filter_method, "/")
dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)

##############################################################################

### Load results

results_padj_list <- list()
truth_list <- list()


for(i in 1:length(simulation_list)){
  
  for(j in 1:length(count_method_list)){
    # i = 1; j = 1
    
    simulation <- simulation_list[i]
    count_method <- count_method_list[j]
    
    comparison_out_tmp <- paste0(simulation, "/", comparison_out)
    out_dir <- paste0(comparison_out_tmp, "/", filter_method, "/", count_method, "_")
    
    if(!file.exists(paste0(out_dir, "cobradata.Rdata")))
      next
    
    load(paste0(out_dir, "cobradata.Rdata"))
    
    message(paste0("Loaded ", out_dir))
    
    results_padj <- cobradata@padj
    rownames(results_padj) <- paste0(simulation, "_", count_method, ".", rownames(cobradata@padj))
    
    truth <-  cobradata@truth
    rownames(truth) <- paste0(simulation, "_", count_method, ".", rownames(cobradata@truth))
    truth$simulation <- simulation
    truth$count_method <- count_method
    
    results_padj_list[[paste0(simulation, "_", count_method)]] <- results_padj
    truth_list[[paste0(simulation, "_", count_method)]] <- truth
    
    
  }
}




results_padj <- rbind.fill(results_padj_list)
rownames(results_padj) <- unlist(lapply(results_padj_list, rownames))

truth <- rbind.fill(truth_list)
rownames(truth) <- unlist(lapply(truth_list, rownames))




table(truth$simulation)
table(truth$count_method)


truth$simulation <- factor(truth$simulation, levels = simulation_list)
levels(truth$simulation)

truth$count_method <- factor(truth$count_method, levels = count_method_list)
levels(truth$count_method)



##############################################################################
### Order for split - simulations per row
##############################################################################


truth$split <- factor(interaction(truth$simulation, truth$count_method), levels = paste(rep(levels(truth$simulation), each = nlevels(truth$count_method)), levels(truth$count_method), sep = ".")
)
levels(truth$split)
table(truth$split)



cobradata <- COBRAData(padj = results_padj, truth = truth)


### FDR TPR stratified

splv <- "split"

cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", splv = splv, aspects = c("fdrtpr", "fdrtprcurve"), onlyshared = FALSE, maxsplit = Inf)

cobraplot <- prepare_data_for_plot(cobraperf, incloverall = FALSE, colorscheme = colors[basemethods(cobraperf)])


colors_df <- colors_df[colors_df$methods %in% basemethods(cobraperf), , drop = FALSE]

cobraplot <- iCOBRA:::reorder_levels(cobraplot, colors_df$methods)



facet_nrow <- length(simulation_list)


### Plot points only

levels(cobraplot@fdrtpr$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtpr$splitval))

xaxisrange <- c(0, 0.7)
yaxisrange <- c(0.4, 0.9)

ggp <- plot_fdrtprcurve(cobraplot, plottype = c("points"), pointsize = 3)
ggp <- ggp + 
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16, color = "darkgrey"), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), strip.text = element_text(size = 11), strip.background = element_rect(colour = "black", fill="white")) + 
  guides(colour = guide_legend(nrow = legend_nrow)) + 
  facet_wrap(~splitval, nrow = facet_nrow) +
  coord_cartesian(xlim = xaxisrange, ylim = yaxisrange)


pdf(paste0(out_dir_plots, "fdrtpr_simulation_count_method", name, ".pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()



### Plot points and a curve

levels(cobraplot@fdrnbrcurve$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrnbrcurve$splitval))
levels(cobraplot@fdrtprcurve$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtprcurve$splitval))


xaxisrange <- c(0, 0.7)
yaxisrange <- c(0.4, 0.9)

ggp <- plot_fdrtprcurve(cobraplot, pointsize = 3)
ggp <- ggp + 
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16, color = "darkgrey"), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), strip.text = element_text(size = 11), strip.background = element_rect(colour = "black", fill="white")) + 
  guides(colour = guide_legend(nrow = legend_nrow)) + 
  facet_wrap(~splitval, nrow = facet_nrow) +
  coord_cartesian(xlim = xaxisrange, ylim = yaxisrange)


pdf(paste0(out_dir_plots, "fdrtprcurve_simulation_count_method", name, ".pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()




##############################################################################
### Order for split - simulations per column
##############################################################################

if(all(simulation_list == c('drosophila_node_nonull','hsapiens_node_nonull')) && all(count_method_list == c('kallisto','kallistofiltered5','htseq','htseqprefiltered5'))){
  
  
  truth$split <- factor(interaction(truth$simulation, truth$count_method), levels = paste(rep(levels(simulation_list), each = nlevels(count_method_list)), levels(count_method_list), sep = ".")
  )
  
  
  truth$split <- interaction(truth$simulation, truth$count_method, lex.order = TRUE)
  
  levels(truth$split)
  
  truth$split <- factor(truth$split, levels = levels(truth$split)[c(1, 3, 5, 7, 2, 4, 6, 8)])
  
  levels(truth$split)
  
  table(truth$split)
  
  
  ### Plot
  
  
  cobradata <- COBRAData(padj = results_padj[, c("dexseq", "drimseq_genewise_grid_none")], truth = truth)
  
  cobradata <- COBRAData(padj = results_padj, truth = truth)
  
  
  ### FDR TPR stratified
  
  splv <- "split"
  
  cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", splv = splv, aspects = c("fdrtpr", "fdrtprcurve"), onlyshared = FALSE, maxsplit = Inf)
  
  cobraplot <- prepare_data_for_plot(cobraperf, incloverall = FALSE, colorscheme = colors[basemethods(cobraperf)])
  
  
  colors_df <- colors_df[colors_df$methods %in% basemethods(cobraperf), , drop = FALSE]
  
  cobraplot <- iCOBRA:::reorder_levels(cobraplot, colors_df$methods)
  
  
  levels(cobraplot@fdrtpr$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtpr$splitval))
  
  facet_nrow <- length(simulation_list)
  
  
  ### Plot points only
  
  xaxisrange <- c(0, 0.7)
  yaxisrange <- c(0.4, 0.9)
  
  ggp <- plot_fdrtprcurve(cobraplot, plottype = c("points"), pointsize = 3)
  ggp <- ggp + 
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16, color = "darkgrey"), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), strip.text = element_text(size = 11), strip.background = element_rect(colour = "black", fill="white")) + 
    guides(colour = guide_legend(nrow = legend_nrow)) + 
    facet_wrap(~splitval, nrow = facet_nrow) +
    coord_cartesian(xlim = xaxisrange, ylim = yaxisrange)
  
  pdf(paste0(out_dir_plots, "fdrtpr_simulation_count_method_order2", name, ".pdf"), width = pdf_width, height = pdf_height)
  print(ggp)
  dev.off()
  
  
  
  ### Plot points and a curve
  
  levels(cobraplot@fdrtprcurve$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtprcurve$splitval))
  
  xaxisrange <- c(0, 0.7)
  yaxisrange <- c(0.4, 0.9)
  
  ggp <- plot_fdrtprcurve(cobraplot, pointsize = 3)
  ggp <- ggp + 
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.text = element_text(size = 16, color = "darkgrey"), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), strip.text = element_text(size = 11), strip.background = element_rect(colour = "black", fill="white")) + 
    guides(colour = guide_legend(nrow = legend_nrow)) + 
    facet_wrap(~splitval, nrow = facet_nrow) +
    coord_cartesian(xlim = xaxisrange, ylim = yaxisrange)
  
  
  
  pdf(paste0(out_dir_plots, "fdrtprcurve_simulation_count_method_order2", name, ".pdf"), width = pdf_width, height = pdf_height)
  print(ggp)
  dev.off()
  
}



##############################################################################
### Order for split - simulations per column
##############################################################################



results <- list()
metadata <- list()


for(i in 1:length(simulation_list)){
  
  for(j in 1:length(count_method_list)){
    # i = 1; j = 1
    
    simulation <- simulation_list[i]
    count_method <- count_method_list[j]
    
    message(paste0(simulation, " ", count_method))
    
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
    
    if(!file.exists(paste0(simulation, "/", results_dir, ".txt")))
      next
    
    rt <- read.table(paste0(simulation, "/", results_dir, ".txt"), header = TRUE, as.is = TRUE)
    head(rt)
    
    colnames(rt) <- c("gene_id", "adj_pvalue")
    rt$pvalue <- rt$adj_pvalue
    head(rt)
    
    method_name <- "dexseq"
    results[[paste(simulation, count_method, method_name, sep = "_")]] <- rt
    metadata[[paste(simulation, count_method, method_name, sep = "_")]] <- data.frame(simulation = simulation, count_method = count_method, method_name = method_name, stringsAsFactors = FALSE)
    
    
    ####################### results from DRIMSeq
    
    results_dir <- paste0(simulation, "/",method_out, "/", count_method, "/", filter_method, "/")
    files <- list.files(path = results_dir, pattern = "_results.txt" )
    files
    
    if(length(files) > 0){
      
      for(ii in 1:length(files)){
        # ii = 1
        method_name <- gsub(pattern = "_results.txt", replacement = "", x = files[ii])
        
        rt <- read.table(paste0(results_dir, files[ii]), header = TRUE, as.is = TRUE)
        head(rt)
        
        rt <- rt[,c("gene_id","adj_pvalue", "pvalue")]
        
        results[[paste(simulation, count_method, method_name, sep = "_")]] <- rt
        metadata[[paste(simulation, count_method, method_name, sep = "_")]] <- data.frame(simulation = simulation, count_method = count_method, method_name = method_name, stringsAsFactors = FALSE)
        
        
      }
      
    }
    
  }
  
}


metadata <- rbind.fill(metadata)


metadata$method_name <- factor(metadata$method_name, levels = colors_df$methods)
metadata$method_name <- factor(metadata$method_name)
metadata$simulation <- factor(metadata$simulation, levels = simulation_list)
metadata$count_method <- factor(metadata$count_method, levels = count_method_list)



source(CAT_function_path)


### Create an index that indicates pairs of methods to be compared
metadata$interaction <- interaction(metadata$simulation, metadata$count_method, lex.order = TRUE)

interaction_levels <- levels(metadata$interaction)

metadata$results_order <- 1:nrow(metadata)

reference_method <- "dexseq"

indx <- lapply(1:nlevels(metadata$interaction), function(i){
  # i = 1
  
  metadata_tmp <- subset(metadata, interaction == interaction_levels[i])
  
  indx1 <- metadata_tmp[metadata_tmp$method_name == reference_method, "results_order"]
  indx2 <- metadata_tmp[!metadata_tmp$method_name == reference_method, "results_order"]
  
  data.frame(indx1 = rep(indx1, length(indx2)), indx2 = indx2)
  
})

indx <- rbind.fill(indx)




### Calculate overlaps for each pair of methods
data_CAT <- lapply(1:nrow(indx), function(j){
  # j = 1
  
  calculateCAT(results1 = results[[indx$indx1[j]]], results2 = results[[indx$indx2[j]]], by = 10)

})


### Update metadata for overlaps in data_Overlaps list
metadata_ov <- metadata[indx$indx2, ]
metadata_ov$method_name <- factor(metadata_ov$method_name)


ggp <- plotCAT(data_CAT, metadata = metadata_ov, plot_var = "method_name", facet_var = c("count_method", "simulation"), plot_colors = colors[levels(metadata_ov$method_name)], plotx = TRUE, reference_color = colors[reference_method])


ggp <- ggp + 
  coord_cartesian(xlim = c(0, 1500), ylim = c(0, 1)) +
  guides(colour = guide_legend(nrow = legend_nrow))


pdf(paste0(out_dir_plots, "cat_simulation_count_method", name, ".pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()

























