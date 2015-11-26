######################################################
## ----- sim5_drimseq_0_3_1_comparison_combined
## <<sim5_drimseq_0_3_1_comparison_combined.R>>

# BioC 3.1
# Created 16 Nov 2015 

##############################################################################

library(iCOBRA)
library(Hmisc)
library(plyr)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_sim5'

simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull','hsapiens_withde_nonull')
count_method_list=c('kallisto','kallistofiltered5','kallistoprefiltered5','htseq','htseqprefiltered5','htseqprefiltered15')
name=''
legend_nrow=1
pdf_width=15
pdf_height=6

##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(simulation_list)
print(count_method_list)
print(legend_nrow)
print(pdf_width)
print(pdf_height)


##############################################################################

setwd(rwd)
method_out <- "drimseq_0_3_1"


### colors

load(paste0(rwd, "/colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)


##############################################################################

### Load results

results_padj_list <- list()
truth_list <- list()


for(i in 1:length(simulation_list)){
  
  for(j in 1:length(count_method_list)){
    # i = 1; j = 1
    
    simulation <- simulation_list[i]
    count_method <- count_method_list[j]
    
    comparison_out <- paste0(simulation, "/drimseq_0_3_1_comparison")
    out_dir <- paste0(comparison_out, "/", count_method, "_")
    
    if(!file.exists(paste0(out_dir, "cobradata.Rdata")))
      next
    
    load(paste0(out_dir, "cobradata.Rdata"))
    
    
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


truth$split <- factor(interaction(truth$simulation, truth$count_method), levels = paste(rep(levels(truth$simulation), each = nlevels(truth$count_method)), levels(truth$count_method), sep = ".")
)
levels(truth$split)
table(truth$split)




### Plot

out_dir <- paste0("drimseq_0_3_1_comparison/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



cobradata <- COBRAData(padj = results_padj[, c("dexseq", "drimseq_genewise_grid_none")], truth = truth)

cobradata <- COBRAData(padj = results_padj, truth = truth)


### FDR TPR stratified

splv <- "split"

cobraperf <- calculate_performance(cobradata, binary_truth = "ds_status", splv = splv, aspects = "fdrtpr", onlyshared = TRUE, maxsplit = Inf)

cobraplot <- prepare_data_for_plot(cobraperf, incloverall = FALSE, colorscheme = colors[basemethods(cobraperf)])


colors_df <- colors_df[colors_df$methods %in% basemethods(cobraperf), , drop = FALSE]

cobraplot@fdrtpr$method <- factor(cobraplot@fdrtpr$method, levels = colors_df$methods)


levels(cobraplot@fdrtpr$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtpr$splitval))

facet_nrow <- length(simulation_list)


ggp <- plot_fdrtprcurve(cobraplot, plottype = c("points"), pointsize = 3, stripsize = 9, xaxisrange = c(0, 0.6), yaxisrange = c(0.4, 1))
ggp <- ggp + 
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = legend_nrow)) + 
  facet_wrap(~splitval, nrow = facet_nrow)


pdf(paste0(out_dir, "fdrtpr_simulation_count_method", name, ".pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()




















