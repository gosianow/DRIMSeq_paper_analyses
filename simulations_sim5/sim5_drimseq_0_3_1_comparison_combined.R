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

##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)


##############################################################################


setwd(rwd)
method_out <- "drimseq_0_3_1"


### Methods

simulation_list=c('drosophila_node_nonull','hsapiens_node_nonull')
count_method_list=c('kallisto','htseq','htseq_prefiltered15')



### Colors
colors_tmp <- read.table(paste0(rwd, "/colors.txt"), header = TRUE, as.is = TRUE)
colors <- colors_tmp$colors
names(colors) <- colors_tmp$methods



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

levels(cobraplot@fdrtpr$splitval) <- gsub(paste0(cobraplot@splv, ":"), "", levels(cobraplot@fdrtpr$splitval))


ggp <- plot_fdrtprcurve(cobraplot, plottype = c("points"), pointsize = 3, stripsize = 9, xaxisrange = c(0, 0.6), yaxisrange = c(0.4, 1))
ggp <- ggp + 
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 2)) + 
  facet_wrap(~splitval, nrow = 2)


pdf(paste0(out_dir, "fdrtpr_simulation_count_method.pdf"), 9, 7)
print(ggp)
dev.off()




















