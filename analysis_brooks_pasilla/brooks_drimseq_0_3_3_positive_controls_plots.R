######################################################
## ----- brooks_drimseq_0_3_3_positive_controls_plots
## <<brooks_drimseq_0_3_3_positive_controls_plots.R>>

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
# models=c('model_full','model_full_paired','model_full_glm')
# method_out='drimseq_0_3_3'
# comparison_out='drimseq_0_3_3_positive_controls'
# ROC_function_path='/home/gosia/R/drimseq_code/help_functions/dm_plotROCx.R'
# valid_path='5_validation/brooks_validated_genes.txt'
# strip_text_size=16


##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(models)
print(count_methods)

##############################################################################

setwd(rwd)

comparison_out <- paste0(comparison_out, "/")

dir.create(comparison_out, showWarnings = FALSE, recursive = TRUE)

dir.create(paste0(comparison_out, "figures/"), showWarnings = FALSE, recursive = TRUE)


##############################################################################
# metadata
##############################################################################


metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors=F, sep="\t", header=T) 

metadata

##############################################################################
# colors
##############################################################################


load(paste0(rwd, "/", comparison_out, "/colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)


##############################################################################
# validated genes
##############################################################################


valid <- read.table(valid_path, header = TRUE, sep = "\t", as.is = TRUE) 

valid

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


### Create results list using adjusted p-values as p-values


results_padj <- lapply(results_padj, function(x){
  
  colnames(x) <- c("gene_id", "pvalue")
  x$adj_pvalue <- x$pvalue
  
  return(x)
  
})



############################################################################
# "ROC" plots using dm_plotROCx.R 
############################################################################

source(ROC_function_path)


### Create "truth" from the validated genes

status <- data.frame(gene_id = results_padj[[1]]$gene_id, status = FALSE, stringsAsFactors = FALSE)
rownames(status) <- status$gene_id

status$status[status$gene_id %in% valid$gene_id] <- TRUE

table(status$status)



data_ROCx <- calculateROCx(results = results_padj, status = status, P = 1, N = 1)


ggp <- plotROCx(data_ROCx, metadata, plot_var = "method_name", facet_var = c("count_method", "model"), plot_colors = colors[levels(metadata$method_name)], plotx = TRUE)


ggp <- ggp + 
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, nrow(valid))) +
  xlab("Number of genes detected as DS") +
  ylab("Number of validated genes detected as DS") +
  guides(colour = guide_legend(nrow = 1)) +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16), strip.text = element_text(size = strip_text_size))



pdf(paste0(comparison_out, "figures/", "number_ds_genes.pdf"), width = 14, height = 10)
print(ggp)
dev.off()







