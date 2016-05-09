######################################################
## <<kim_drimseq_0_3_3_comparison.R>>

# BioC 3.2
# Created 15 Jan 2015 
# Modified 14 Apr 2016

##############################################################################

Sys.time()

#######################################################

library(ggplot2)
library(iCOBRA)
library(DRIMSeq)
library(limma)


##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/kim_adenocarcinoma'
# count_method=c('htseq','kallisto')[2]
# model=c('model_full','model_full_glm','model_null_normal1','model_null_tumor1')[2]
# method_out='drimseq_0_3_3'
# comparison_out='drimseq_0_3_3_comparison'
# text_size=14
# legend_size=12


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
    
    text_size_disp <- 26
    
    ggp <- plotDispersion(d) +
      coord_cartesian(ylim = log10(c(splineDisp[1], splineDisp[disp_grid_length]))) +
      theme(axis.text = element_text(size = text_size_disp), axis.title = element_text(size = text_size_disp, face = "bold"), legend.text = element_text(size = text_size_disp), legend.title = element_text(size = text_size_disp, face = "bold"))
    
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
    
    ### Plot p-values
    ggp <- plotTest(d) +
      theme(axis.text = element_text(size = text_size_disp), axis.title = element_text(size = text_size_disp, face = "bold"), plot.title = element_text(size = text_size_disp))
    
    pdf(paste0(res_path, gsub("_d.Rdata", "", files[i]), "_hist_pvalues.pdf"))
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

clolors <- colors[keep_methods]
colors_df <- colors_df[keep_methods, , drop = FALSE]


#######################################################
### use iCOBRA
#######################################################

summary <- data.frame(model = model, count_method = count_method, ds_method = colnames(results_padj))

if(nrow(summary) == 1){
  
  summary$counts_genes_all <- summary$counts_genes_all_overlap <- sum(!is.na(results_padj))
  summary$counts_genes_ds <- summary$counts_genes_ds_overlap <- sum(results_padj[!is.na(results_padj), ] < 0.05)
  
}else{
  
  stopifnot("dexseq" %in% colnames(results_padj))
  
  cobradata <- COBRAData(padj = results_padj)
  
  ### Count how many genes are tested in total
  
  cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 1.1)
  
  overlap <- cobraperf@overlap 
  overlap[is.na(overlap)] <- 0
  
  summary$counts_genes_all <- apply(overlap, 2, sum)
  
  ref_dexseq <- overlap[, "dexseq"]
  
  summary$counts_genes_all_overlap <- apply(overlap, 2, function(i){
    sum((ref_dexseq + i) == 2)
  })
  
  
  ## Plot overlaps
  for(i in 2:length(basemethods(cobraperf))){
    # i = 2
    
    cobraplot <- prepare_data_for_plot(cobraperf, keepmethods = c("dexseq", basemethods(cobraperf)[i]), colorscheme = c(colors[c("dexseq", basemethods(cobraperf)[i])]), incltruth = FALSE)
    
    pdf(paste0(out_dir, "/venn_", basemethods(cobraperf)[i], "_all.pdf"))
    plot_overlap(cobraplot, cex=c(1.2,1,0.7))
    dev.off()
    
    pdf(paste0(out_dir, "/upset_", basemethods(cobraperf)[i], "_all.pdf"), width = 10)
    plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 14, sets = c("dexseq", basemethods(cobraperf)[i]), sets.bar.color = colors[c("dexseq", basemethods(cobraperf)[i])])
    dev.off()
    
  }
  
  
  ### Count how many genes are siginificant
  
  cobraperf <- calculate_performance(cobradata, aspects = "overlap", thr_venn = 0.05)
  
  
  overlap <- cobraperf@overlap 
  overlap[is.na(overlap)] <- 0
  
  summary$counts_genes_ds <- apply(overlap, 2, sum)
  
  ref_dexseq <- overlap[, "dexseq"]
  
  summary$counts_genes_ds_overlap <- apply(overlap, 2, function(i){
    sum((ref_dexseq + i) == 2)
  })
  
  
  ## Plot overlaps
  for(i in 2:length(basemethods(cobraperf))){
    # i = 2
    
    cobraplot <- prepare_data_for_plot(cobraperf, keepmethods = c("dexseq", basemethods(cobraperf)[i]), colorscheme = c(colors[c("dexseq", basemethods(cobraperf)[i])]), incltruth = FALSE)
    
    pdf(paste0(out_dir, "/venn_", basemethods(cobraperf)[i], ".pdf"))
    plot_overlap(cobraplot, cex=c(1.2,1,0.7))
    dev.off()
    
    pdf(paste0(out_dir, "/upset_", basemethods(cobraperf)[i], ".pdf"), width = 10)
    plot_upset(cobraplot, order.by = "degree", empty.intersections = "on", name.size = 14, sets = c("dexseq", basemethods(cobraperf)[i]), sets.bar.color = colors[c("dexseq", basemethods(cobraperf)[i])])
    dev.off()
    
  }
  
  if(model == "model_full"){
    
    ### load metadata
    metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors = FALSE, sep="\t", header=TRUE) 
    
    metadata_org <- metadata
    
    count_dir <- paste0("2_counts/", count_method, "/")
    
    
    ### load counts
    counts_list <- lapply(1:length(metadata_org$sampleName), function(i){
      # i = 1
      cts <- read.table(paste0(count_dir, metadata_org$sampleName[i], ".counts"), header = FALSE, as.is = TRUE)
      colnames(cts) <- c("group_id", metadata_org$sampleName[i])  
      return(cts)
    })
    
    counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)
    counts <- counts[!grepl(pattern = "_", counts$group_id),]
    
    
    ### Prepare data
    group_split <- strsplit2(counts[,1], ":")
    counts <- counts[, -1]
    ### order the samples like in metadata!!!
    counts <- counts[, metadata_org$sampleName]
    
    
    ### Plot distributions of differenet gene characteristics for unique genes
    bmethds <- basemethods(cobraperf)
    
    
    ### Mean expression
    
    mean_expression <- by(counts, factor(group_split[, 1]), function(x){
      mean(colSums(x, na.rm = TRUE), na.rm = TRUE)
    }, simplify = FALSE)
    
    mean_expression <- unlist(mean_expression)
    
    
    ### Number of expressed features per gene 
    
    nbr_features<- by(counts, factor(group_split[, 1]), function(x){
      x <- as.matrix(x)
      sum(rowSums(x > 10, na.rm = TRUE) > 2, na.rm = TRUE)
    }, simplify = FALSE)
    
    nbr_features <- unlist(nbr_features)
    
    
    
    for(i in grep("drimseq", bmethds)){
      # i = 2
      
      all_genes <- rownames(results_padj)
      
      m1 <- "dexseq"
      m2 <- bmethds[i]
      
      genes_sign_m1 <- rownames(results_padj[results_padj[, m1] < 0.05 & !is.na(results_padj[, m1]), , drop = FALSE])
      genes_sign_m2 <- rownames(results_padj[results_padj[, m2] < 0.05 & !is.na(results_padj[, m2]), ,drop = FALSE])
      
      genes_sign_overlap <- intersect(genes_sign_m1, genes_sign_m2)
      genes_sign_m1_unique <- setdiff(genes_sign_m1, genes_sign_overlap)
      genes_sign_m2_unique <- setdiff(genes_sign_m2, genes_sign_overlap)
      
      
      ggdf <- data.frame(mean_expression = c(mean_expression[all_genes], mean_expression[genes_sign_m1], mean_expression[genes_sign_m2]), 
        group = c(rep("all_genes", length(all_genes)), rep(m1, length(genes_sign_m1)), rep(m2, length(genes_sign_m2))), stringsAsFactors = FALSE)
      
      ggdf <- ggdf[complete.cases(ggdf), ]
      
      ggdf$group <- factor(ggdf$group, levels = c("all_genes", m1, m2), labels = paste0(c("all_genes", m1, m2), " (", c(length(all_genes), length(genes_sign_m1), length(genes_sign_m2)), ")"))
      
      ggdf <- ggdf[ggdf$mean_expression > 0, ,drop = FALSE]
      
      
      ggp <- ggplot(ggdf, aes(x = log10(mean_expression), color = group, group = group)) +
        geom_density(size = 2) +
        theme_bw() +
        xlab("Log10 of gene mean expression") +
        theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = legend_size), legend.title = element_blank(), legend.position = "bottom") +
        scale_color_manual(values = as.character(c("grey", colors[c(m1, m2)]))) +
        guides(color = guide_legend(nrow = 2))
      
      pdf(paste0(out_dir, "characteristics_mean_expr_", bmethds[i], ".pdf"))
      print(ggp)
      dev.off()
      
      
      ### Unique genes
      ggdf <- data.frame(mean_expression = c(mean_expression[all_genes], mean_expression[genes_sign_overlap], mean_expression[genes_sign_m1_unique], mean_expression[genes_sign_m2_unique]), 
        group = c(rep("all_genes", length(all_genes)), rep("overlap", length(genes_sign_overlap)), rep(m1, length(genes_sign_m1_unique)), rep(m2, length(genes_sign_m2_unique))), stringsAsFactors = FALSE)
      
      ggdf <- ggdf[complete.cases(ggdf), ]
      
      ggdf$group <- factor(ggdf$group, levels = c("all_genes", "overlap", m1, m2), labels = paste0(c("all_genes", "overlap", m1, m2), " (", c(length(all_genes), length(genes_sign_overlap), length(genes_sign_m1_unique), length(genes_sign_m2_unique)), ")"))
      
      ggdf <- ggdf[ggdf$mean_expression > 0, ,drop = FALSE]
      
      
      ggp <- ggplot(ggdf, aes(x = log10(mean_expression), color = group, group = group)) +
        geom_density(size = 2) +
        theme_bw() +
        xlab("Log10 of gene mean expression") +
        theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = legend_size), legend.title = element_blank(), legend.position = "bottom") +
        scale_color_manual(values = as.character(c("grey", "gray50", colors[c(m1, m2)]))) +
        guides(color = guide_legend(nrow = 2))
      
      pdf(paste0(out_dir, "characteristics_mean_expr_", bmethds[i], "2.pdf"))
      print(ggp)
      dev.off()
      
      
      
      
      ggdf <- data.frame(nbr_features = nbr_features[c(all_genes, genes_sign_m1, genes_sign_m2)], 
        group = rep(c("all_genes", m1, m2), c(length(all_genes), length(genes_sign_m1), length(genes_sign_m2))),
        stringsAsFactors = FALSE)
      
      ggdf <- ggdf[complete.cases(ggdf), ]
      
      ggdf$group <- factor(ggdf$group, levels = c("all_genes", m1, m2), labels = paste0(c("all_genes", m1, m2), " (", c(length(all_genes), length(genes_sign_m1), length(genes_sign_m2)), ")"))
      
      
      ggp <- ggplot(ggdf, aes(x = nbr_features, color = group, group = group)) +
        geom_density(size = 2) +
        theme_bw() +
        xlab("Number of expressed features") +
        theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = legend_size), legend.title = element_blank(), legend.position = "bottom") +
        scale_color_manual(values = as.character(c("grey", colors[c(m1, m2)])))+
        guides(color = guide_legend(nrow = 2))
      
      pdf(paste0(out_dir, "characteristics_nbr_features_", bmethds[i], ".pdf"))
      print(ggp)
      dev.off()
      
      
      ### Unique genes
      ggdf <- data.frame(nbr_features = nbr_features[c(all_genes, genes_sign_overlap, genes_sign_m1_unique, genes_sign_m2_unique)], 
        group = rep(c("all_genes", "overlap", m1, m2), c(length(all_genes), length(genes_sign_overlap), length(genes_sign_m1_unique), length(genes_sign_m2_unique))),
        stringsAsFactors = FALSE)
      
      ggdf <- ggdf[complete.cases(ggdf), ]
      
      ggdf$group <- factor(ggdf$group, levels = c("all_genes", "overlap", m1, m2), labels = paste0(c("all_genes", "overlap", m1, m2), " (", c(length(all_genes), length(genes_sign_overlap), length(genes_sign_m1_unique), length(genes_sign_m2_unique)), ")"))
      
      
      ggp <- ggplot(ggdf, aes(x = nbr_features, color = group, group = group)) +
        geom_density(size = 2) +
        theme_bw() +
        xlab("Number of expressed features") +
        theme(axis.text = element_text(size = text_size), axis.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = legend_size), legend.title = element_blank(), legend.position = "bottom") +
        scale_color_manual(values = as.character(c("grey", "gray50", colors[c(m1, m2)])))+
        guides(color = guide_legend(nrow = 2))
      
      pdf(paste0(out_dir, "characteristics_nbr_features_", bmethds[i], "2.pdf"))
      print(ggp)
      dev.off()
      
      
    }
    
    
  }
  
  
}

summary

write.table(summary, file = paste0(comparison_out, model, "_", count_method, "_summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




sessionInfo()
