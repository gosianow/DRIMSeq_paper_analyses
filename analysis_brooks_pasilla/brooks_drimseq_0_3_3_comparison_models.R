##############################################################################
## <<brooks_drimseq_0_3_3_comparison_models.R>>

# BioC 3.2
# Created 14 Apr 2016 

##############################################################################
Sys.time()
##############################################################################

library(ggplot2)
library(iCOBRA)
library(DRIMSeq)


##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/brooks_pasilla'
# count_method='kallisto'
# ds_method='drimseq_genewise_grid_trended'
# model_list=c('model_full','model_full_paired')
# method_out='drimseq_0_3_3'
# comparison_out='drimseq_0_3_3_comparison'

##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)


##############################################################################

setwd(rwd)

comparison_out <- paste0(comparison_out, "/")

out_dir <- paste0(comparison_out,  "between_models", "/", count_method, "/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

### colors

load(paste0(rwd, "/", comparison_out, "colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)


#######################################################
# merge results for iCOBRA
#######################################################


results_padj <- list()


####################### results from DEXSeq

if(grepl("dexseq", ds_method)){
  
  for(model in model_list){
    
    rt <- read.table(paste0("4_results/dexseq_1_10_8/", model,"/", count_method, "/dexseq_gene_results.txt"), header = TRUE, as.is = TRUE)
    head(rt)
    
    colnames(rt) <- c("gene_id", model)
    head(rt)
    
    results_padj[[model]] <- rt
    
    
  }
  
}



####################### results from DRIMSeq


if(grepl("drimseq", ds_method)){
  
  for(model in model_list){
    
    res_path <- paste0(method_out, "/",  model, "/", count_method, "/")
    files <- list.files(path = res_path, pattern = paste0(ds_method, "_results.txt"))
    files
    
    if(length(files) > 0){
      for(i in 1:length(files)){
        # i = 1
        method_name <- gsub(pattern = "_results.txt", replacement = "", x = files[i])
        
        rt <- read.table(paste0(res_path, files[i]), header = TRUE, as.is = TRUE)
        head(rt)
        
        rt <- rt[,c("gene_id","adj_pvalue")]
        colnames(rt) <- c("gene_id", model)
        
        results_padj[[model]] <- rt 
        
      }
    }
    
  }
}


####################### merge results


results_padj <- Reduce(function(...) merge(..., by = "gene_id", all=TRUE, sort = FALSE), results_padj)
rownames(results_padj) <- results_padj$gene_id

colnames(results_padj)


results_padj <- results_padj[, !grepl("gene_id", colnames(results_padj)), drop = FALSE]




#######################################################
### use iCOBRA
#######################################################

### standard ggplot scheme

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}


cobradata <- COBRAData(padj = results_padj)

### Count how many genes are tested in total

cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 1.1)

bmethods <- basemethods(cobraperf)


cobraplot <- prepare_data_for_plot(cobraperf, incltruth = FALSE)

pdf(paste0(out_dir, "/venn_", ds_method, "_all.pdf"))
plot_overlap(cobraplot, cex=c(1.2,1,0.7))
dev.off()

pdf(paste0(out_dir, "/upset_", ds_method, "_all.pdf"))
plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 14, sets = bmethods, sets.bar.color = gg_color_hue(length(bmethods)))
dev.off()



### Count how many genes are siginificant
cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 0.05)

cobraplot <- prepare_data_for_plot(cobraperf, incltruth = FALSE)

pdf(paste0(out_dir, "/venn_", ds_method, ".pdf"))
plot_overlap(cobraplot, cex=c(1.2,1,0.7))
dev.off()

pdf(paste0(out_dir, "/upset_", ds_method, ".pdf"))
plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 14, sets = bmethods, sets.bar.color = gg_color_hue(length(bmethods)))
dev.off()



















































