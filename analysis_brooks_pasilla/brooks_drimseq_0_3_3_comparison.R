##############################################################################
## <<brooks_drimseq_0_3_3_comparison_run.R>>

# BioC 3.2
# Created 15 Jan 2015 

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
# model='model_full'
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

print(rwd)
print(model)
print(count_method)


##############################################################################

setwd(rwd)

comparison_out <- paste0(comparison_out, "/")

out_dir <- paste0(comparison_out,  model, "/", count_method, "/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

### colors

load(paste0(rwd, "/", comparison_out, "colors.Rdata"))
colors
colors_df

colors_df$methods <- as.character(colors_df$methods)


#######################################################
# Redo the dispersion plots
#######################################################


if(model != "model_full_glm"){
  
  res_path <- paste0(method_out, "/",  model, "/", count_method, "/")
  files <- list.files(path = res_path, pattern = "_d.Rdata" )
  files
  
  
  ### For fitting smooth line, do not take into account the boudry grid points
  disp_init <- common_disp <- as.numeric(read.table(paste0(res_path, "common_dispersion.txt")))
  
  disp_grid_length = 21 
  disp_grid_range = c(-10, 10)
  splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length)
  splineDisp <- disp_init * 2^splinePts
  
  
  for(i in grep("_grid_", files)){
    # i = 3
    
    load(paste0(res_path, files[i]))
    
    ggp <- plotDispersion(d) +
      coord_cartesian(ylim = log10(c(splineDisp[1], splineDisp[disp_grid_length])))
    
    if(grepl("_grid_none_", files[i])){
      
      ggdf <- ggp$data
      not_boundry <- ggdf$dispersion < log10(splineDisp[disp_grid_length - 3]) & ggdf$dispersion > log10(splineDisp[1])
      ggdf <- ggdf[not_boundry, ]
      
      ggp <- ggp + 
        geom_smooth(data = ggdf, color = "black", span = 0.1)
      
    }
    
    pdf(paste0(res_path, gsub("_d.Rdata", "", files[i]), "_dispersion_vs_mean.pdf"))
    print(ggp)
    dev.off()
    
    
  }
  
  
}



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
files

if(length(files) > 0){
  for(i in 1:length(files)){
    # i = 1
    method_name <- gsub(pattern = "_results.txt", replacement = "", x = files[i])
    
    rt <- read.table(paste0(res_path, files[i]), header = TRUE, as.is = TRUE)
    head(rt)
    
    rt <- rt[,c("gene_id","adj_pvalue")]
    colnames(rt) <- c("gene_id", method_name)
    
    results_padj[[method_name]] <- rt 
    
  }
}

####################### Merge results

results_padj <- Reduce(function(...) merge(..., by = "gene_id", all=TRUE, sort = FALSE), results_padj)
rownames(results_padj) <- results_padj$gene_id


results_padj <- results_padj[, colnames(results_padj) %in% colors_df$methods, drop = FALSE]

keep_methods <- colors_df$methods %in% colnames(results_padj)

colors <- colors[keep_methods]
colors_df <- colors_df[keep_methods, , drop = FALSE]


#######################################################
### use iCOBRA
#######################################################

summary <- data.frame(model = model, count_method = count_method, ds_method = colnames(results_padj))
summary

if(nrow(summary) == 1){
  
  summary$counts_genes_all <- summary$counts_genes_all_overlap <- sum(!is.na(results_padj))
  summary$counts_genes_ds <- summary$counts_genes_ds_overlap <- sum(results_padj[!is.na(results_padj), ] < 0.05)
  
}else{
  
  stopifnot("dexseq" %in% colnames(results_padj))
  
  cobradata <- COBRAData(padj = results_padj)
  
  ### Count how many genes are tested in total
  
  cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 1.1)
  
  overlap <- !is.na(cobraperf@overlap)
  
  summary$counts_genes_all <- apply(overlap, 2, sum)
  
  ref_method <- overlap[, "dexseq"]
  
  summary$counts_genes_all_overlap <- apply(overlap, 2, function(i){
    sum((ref_method + i) == 2)
  })
  
  
  ## Plot overlaps
  for(i in grep("drimseq", basemethods(cobraperf))){
    # i = 2
    
    cobraplot <- prepare_data_for_plot(cobraperf, keepmethods = c("dexseq", basemethods(cobraperf)[i]), colorscheme = c(colors[c("dexseq", basemethods(cobraperf)[i])]), incltruth = FALSE)
    
    pdf(paste0(out_dir, "/venn_", basemethods(cobraperf)[i], "_all.pdf"))
    plot_overlap(cobraplot, cex=c(1.2,1,0.7))
    dev.off()
    
    pdf(paste0(out_dir, "/upset_", basemethods(cobraperf)[i], "_all.pdf"))
    plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 10)
    dev.off()
    
  }
  
  
  ### Count how many genes are siginificant
  cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 0.05)
  
  
  overlap <- cobraperf@overlap 
  overlap[is.na(overlap)] <- 0
  
  summary$counts_genes_ds <- apply(overlap, 2, sum)
  
  ref_method <- overlap[, "dexseq"]
  
  summary$counts_genes_ds_overlap <- apply(overlap, 2, function(i){
    sum((ref_method + i) == 2)
  })
  
  ## Plot overlaps
  for(i in grep("drimseq", basemethods(cobraperf))){
    # i = 2
    
    cobraplot <- prepare_data_for_plot(cobraperf, keepmethods = c("dexseq", basemethods(cobraperf)[i]), colorscheme = c(colors[c("dexseq", basemethods(cobraperf)[i])]), incltruth = FALSE)
    
    pdf(paste0(out_dir, "/venn_", basemethods(cobraperf)[i], ".pdf"))
    plot_overlap(cobraplot, cex=c(1.2,1,0.7))
    dev.off()
    
    pdf(paste0(out_dir, "/upset_", basemethods(cobraperf)[i], ".pdf"))
    plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 10)
    dev.off()
    
  }
  
  
}


summary

write.table(summary, file = paste0(comparison_out, model, "_", count_method, "_summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
















































