######################################################
## ----- brooks_drimseq_0_3_3_comparison_plots
## <<brooks_drimseq_0_3_3_comparison_plots.R>>

# BioC 3.2
# Created 24 Feb 2016 

##############################################################################
Sys.time()
##############################################################################

library(plyr)
library(ggplot2)
library(rtracklayer)
library(DRIMSeq)
library(iCOBRA)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/brooks_pasilla'
# count_methods=c('kallisto','kallistofiltered5','htseq','htseqprefiltered5')
# models=c('model_full','model_full_paired')
# method_out='drimseq_0_3_3'
# comparison_out='drimseq_0_3_3_comparison'
# CAT_function_path='/home/gosia/R/drimseq_code/help_functions/dm_plotCAT.R'
# text_size=18
# legend_size=16
# strip_size=16


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

dir.create(comparison_out, showWarnings = FALSE, recursive = TRUE)

dir.create(paste0(comparison_out, "figures/"), showWarnings = FALSE, recursive = TRUE)


##############################################################################
# colors
##############################################################################


load(paste0(rwd, "/", comparison_out, "/colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)



#######################################################
# merge results 
#######################################################


results_padj <- list()
metadata <- list()



for(model in models){
  
  for(count_method in count_methods){
    
    
    ####################### results from DEXSeq
    
    rt <- read.table(paste0("4_results/dexseq_1_10_8/", model,"/", count_method, "/dexseq_gene_results.txt"), header = TRUE, as.is = TRUE)
    head(rt)
    
    colnames(rt) <- c("gene_id", "dexseq")
    head(rt)
    
    method_name <- "dexseq"
    results_padj[[paste(model, count_method, method_name, sep = "_")]] <- rt
    metadata[[paste(model, count_method, method_name, sep = "_")]] <- data.frame(model = model, count_method = count_method, method_name = method_name, stringsAsFactors = FALSE)
    

    ####################### results from DRIMSeq
    
    res_path <- paste0(method_out, "/",  model, "/", count_method, "/")
    files <- list.files(path = res_path, pattern = "_results.txt" )
    files
    
    if(length(files) > 0){
      for(i in 1:length(files)){
        # i = 1
        method_name <- gsub(pattern = "_results.txt", replacement = "", x = files[i])
        
        rt <- read.table(paste0(res_path, files[i]), header = TRUE, as.is = TRUE)
        head(rt)
        
        rt <- rt[,c("gene_id","adj_pvalue")]
        colnames(rt) <- c("gene_id", method_name)
        
        results_padj[[paste(model, count_method, method_name, sep = "_")]] <- rt 
        metadata[[paste(model, count_method, method_name, sep = "_")]] <- data.frame(model = model, count_method = count_method, method_name = method_name, stringsAsFactors = FALSE)
        
      }
    }
    
  }
  
}


metadata <- rbind.fill(metadata)


metadata$method_name <- factor(metadata$method_name, levels = colors_df$methods)
metadata$method_name <- factor(metadata$method_name)
metadata$model <- factor(metadata$model, levels = models)
metadata$count_method <- factor(metadata$count_method, levels = count_methods)

### Data frames in results_padj must have gene_id, pvalue and adj_pvalue columns
### Create results list using adjusted p-values as p-values


results_padj <- lapply(results_padj, function(x){
  
  colnames(x) <- c("gene_id", "pvalue")
  x$adj_pvalue <- x$pvalue
  
  return(x)
  
})



############################################################################
# CAT plots
############################################################################


source(CAT_function_path)


metadata$interaction <- interaction(metadata$model, metadata$count_method, lex.order = TRUE)

interaction_levels <- levels(metadata$interaction)

metadata$results_order <- 1:nrow(metadata)

reference_method <- "dexseq"

### Index that indicates pairs of methods to be compared
indx <- lapply(1:nlevels(metadata$interaction), function(i){
  # i = 1
  
  metadata_tmp <- subset(metadata, interaction == interaction_levels[i])
  
  indx1 <- metadata_tmp[metadata_tmp$method_name == reference_method, "results_order"]
  indx2 <- metadata_tmp[!metadata_tmp$method_name == reference_method, "results_order"]
  
  data.frame(indx1 = rep(indx1, length(indx2)), indx2 = indx2)
  
})

indx <- rbind.fill(indx)


data_CAT <- lapply(1:nrow(indx), function(j){
  # j = 1
  
  calculateCAT(results1 = results_padj[[indx$indx1[j]]], results2 = results_padj[[indx$indx2[j]]], by = 10)
  
  
})


### Metadata for overlaps in data_Overlaps list
metadata_ov <- metadata[indx$indx2, ]
metadata_ov$method_name <- factor(metadata_ov$method_name)


ggp <- plotCAT(data_CAT, metadata = metadata_ov, plot_var = "method_name", facet_var = c("count_method", "model"), plot_colors = colors[levels(metadata_ov$method_name)], plotx = TRUE, reference_color = colors[reference_method])



ggp <- ggp + 
  coord_cartesian(xlim = c(0, 700), ylim = c(0, 1)) +
  theme(legend.position = "bottom", axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = legend_size), strip.text = element_text(size = strip_size), strip.background = element_rect(colour = "black", fill="white"))


pdf(paste0(comparison_out, "figures/", "cat.pdf"), width = 14, height = 7)
print(ggp)
dev.off()




