######################################################
## ----- brooks_drimseq_0_3_3_run
## <<brooks_drimseq_0_3_3_run.R>>

# BioC 3.2
# Created 29 Dec 2015 

##############################################################################
Sys.time()
##############################################################################

library(BiocParallel)
library(DRIMSeq)
library(limma)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/brooks_pasilla/'
# workers=4
# count_method=c('htseq','kallisto')[2]
# model=c('model_full','model_full_paired','model_null1','model_null2','model_null3')[1]
# dispersion_common=TRUE
# results_common=TRUE
# disp_mode_list=c('grid','grid','optimize','optim','constrOptim')
# disp_moderation_list=c('none','common','none','none','none')
# disp_prior_df=0.1

##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(workers)
print(count_method)
print(model)
print(dispersion_common)
print(results_common)
print(disp_mode_list)
print(disp_moderation_list)


##############################################################################

setwd(rwd)
method_out <- "drimseq_0_3_3"


if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}


##########################################################################
# load metadata
##########################################################################

metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors = FALSE, sep="\t", header=TRUE) 

metadata_org <- metadata


##########################################################################
# DRIMSeq analysis
##########################################################################

out_dir <- paste0(method_out, "/",  model, "/", count_method, "/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


count_dir <- paste0("2_counts/", count_method, "/")


### load counts
counts_list <- lapply(1:length(metadata_org$sampleName), function(i){
  # i = 1
  cts <- read.table(paste0(count_dir, metadata_org$sampleName[i], ".txt"), header = FALSE, as.is = TRUE)
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


switch(
  model, 
  
  model_full = {
    
    metadata <- metadata_org
    all(colnames(counts) == metadata$sampleName)
    
    d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sampleName, group = metadata$condition)
    save(d, file = paste0(out_dir, "d.Rdata"))
    ### Filtering
    table(samples(d)$group)
    d <- dmFilter(d, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0, max_features = Inf)

  },
  
  model_full_paired = {
    
    counts <- counts[, metadata_org$LibraryLayout == "PAIRED"]
    metadata <- metadata_org[metadata_org$LibraryLayout == "PAIRED", ]
    all(colnames(counts) == metadata$sampleName)
    
    d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sampleName, group = metadata$condition)
    save(d, file = paste0(out_dir, "d.Rdata"))
    ### Filtering
    table(samples(d)$group)
    d <- dmFilter(d, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0, max_features = Inf)

  },
  
  model_null1 = {

    counts <- counts[, metadata_org$condition == "CTL"]
    metadata <- metadata_org[metadata_org$condition == "CTL", ]
    all(colnames(counts) == metadata$sampleName)
    
    d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sampleName, group = rep(c("C1", "C2"), 2))
    save(d, file = paste0(out_dir, "d.Rdata"))
    ### Filtering
    table(samples(d)$group)
    d <- dmFilter(d, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0, max_features = Inf)

  },
  
  model_null2 = {

    counts <- counts[, metadata_org$condition == "CTL"]
    metadata <- metadata_org[metadata_org$condition == "CTL", ]
    all(colnames(counts) == metadata$sampleName)
    
    d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sampleName, group = rep(c("C1", "C2"), each = 2))
    save(d, file = paste0(out_dir, "d.Rdata"))
    ### Filtering
    table(samples(d)$group)
    d <- dmFilter(d, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0, max_features = Inf)

  },
  
  model_null3 = {
    
    counts <- counts[, metadata_org$condition == "CTL"]
    metadata <- metadata_org[metadata_org$condition == "CTL", ]
    all(colnames(counts) == metadata$sampleName)
    
    d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sampleName, group = c("C1", "C2", "C2", "C1"))
    save(d, file = paste0(out_dir, "d.Rdata"))
    ### Filtering
    table(samples(d)$group)
    d <- dmFilter(d, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0, max_features = Inf)
    
  }
)

plotData(d, out_dir = out_dir)



### DRIMSeq pipelines : common_dispersion
if(dispersion_common){
  
  disp <- "common"
  out_name <- paste0(out_dir, "/drimseq_", disp, "_")
  
  d <- dmDispersion(d, mean_expression = TRUE, common_dispersion = TRUE, genewise_dispersion = FALSE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 0.1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)
  
  common_disp <- common_dispersion(d)
  common_disp
  write.table(common_disp, paste0(out_dir, "common_dispersion.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  
  if(results_common){
    
    d <- dmFit(d, dispersion = "common_dispersion", BPPARAM = BPPARAM)
    
    d <- dmTest(d, BPPARAM = BPPARAM)
    
    plotTest(d, out_dir = out_name)
    
    res <- results(d)
    
    write.table(res, paste0(out_name, "results.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    save(d, file = paste0(out_name, "d.Rdata"))
    
  }
  
}



### DRIMSeq pipelines : genewise_dispersion

common_disp <- as.numeric(read.table(paste0(out_dir, "common_dispersion.txt")))
common_disp


for(i in 1:length(disp_mode_list)){
  # i = 3
  print(i)
  
  load(paste0(out_dir, "d.Rdata"))
  
  disp <- "genewise"
  disp_mode <- disp_mode_list[i]
  disp_moderation <- disp_moderation_list[i]
  
  if(disp_mode == "grid")
    out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_", disp_moderation, "_")
  else
    out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_")
  
  # genewise dispersion
  d <- dmDispersion(d, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_mode = disp_mode, disp_init = common_disp, disp_moderation = disp_moderation, disp_prior_df = disp_prior_df, verbose = TRUE, BPPARAM = BPPARAM)
  common_dispersion(d) <- common_disp
  
  plotDispersion(d, out_dir = out_name)
  
  d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
  
  d <- dmTest(d, BPPARAM = BPPARAM)
  
  plotTest(d, out_dir = out_name)
  
  res <- results(d)
  
  write.table(res, paste0(out_name, "results.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  save(d, file = paste0(out_name, "d.Rdata"))
  
  rm("d")
  gc()
  
  
}


sessionInfo()







