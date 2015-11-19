######################################################
## ----- brooks_drimseq_0_3_1_run
## <<brooks_drimseq_0_3_1_run.R>>

# BioC 3.1
# Created 16 Nov 2015 

##############################################################################

library(DRIMSeq)
library(limma)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/brooks_pasilla/'
workers=4
count_method=c('htseq','kallisto')[2]
model=c('model_full','model_full_paired','model_null1','model_null2','model_null3')[1]
dispersion_common=TRUE
results_common=TRUE
disp_mode_list=c('grid','grid','optimize','optim','constrOptim')
disp_moderation_list=c('none','common','none','none','none')


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
method_out <- "drimseq_0_3_1"

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



### Prepare data
counts <- read.table(paste0("2_counts/", count_method, "/", count_method, "_counts.txt"), header = TRUE, as.is = TRUE)
counts <- counts[!grepl(pattern = "_", counts$group_id),]
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
    
    ### Filtering
    table(samples(d)$group)
    d <- dmFilter(d, min_samps_gene_expr = 3, min_samps_feature_prop = 3, min_feature_prop = 0.01)

  },
  
  model_full_paired = {
    
    counts <- counts[, metadata_org$LibraryLayout == "PAIRED"]
    metadata <- metadata_org[metadata_org$LibraryLayout == "PAIRED", ]
    all(colnames(counts) == metadata$sampleName)
    
    d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sampleName, group = metadata$condition)
    
    ### Filtering
    table(samples(d)$group)
    d <- dmFilter(d, min_samps_gene_expr = 2, min_samps_feature_prop = 2, min_feature_prop = 0.01)

  },
  
  model_null1 = {

    counts <- counts[, metadata_org$condition == "CTL"]
    metadata <- metadata_org[metadata_org$condition == "CTL", ]
    all(colnames(counts) == metadata$sampleName)
    
    d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sampleName, group = rep(c("C1", "C2"), 2))
    
    ### Filtering
    table(samples(d)$group)
    d <- dmFilter(d, min_samps_gene_expr = 2, min_samps_feature_prop = 2, min_feature_prop = 0.01)

  },
  
  model_null2 = {

    counts <- counts[, metadata_org$condition == "CTL"]
    metadata <- metadata_org[metadata_org$condition == "CTL", ]
    all(colnames(counts) == metadata$sampleName)
    
    d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sampleName, group = rep(c("C1", "C2"), each = 2))
    
    ### Filtering
    table(samples(d)$group)
    d <- dmFilter(d, min_samps_gene_expr = 2, min_samps_feature_prop = 2, min_feature_prop = 0.01)

  },
  
  model_null3 = {
    
    counts <- counts[, metadata_org$condition == "CTL"]
    metadata <- metadata_org[metadata_org$condition == "CTL", ]
    all(colnames(counts) == metadata$sampleName)
    
    d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sampleName, group = c("C1", "C2", "C2", "C1"))
    
    ### Filtering
    table(samples(d)$group)
    d <- dmFilter(d, min_samps_gene_expr = 2, min_samps_feature_prop = 2, min_feature_prop = 0.01)
    
  }
)

plotData(d, out_dir = out_dir)
save(d, file = paste0(out_dir, "d.Rdata"))



### DRIMSeq pipelines : common_dispersion
if(dispersion_common){
  
  disp <- "common"
  out_name <- paste0(out_dir, "/drimseq_", disp, "_")
  
  d <- dmDispersion(d, mean_expression = TRUE, common_dispersion = TRUE, genewise_dispersion = FALSE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  common_disp <- common_dispersion(d)
  common_disp
  write.table(common_disp, paste0(out_dir, "common_dispersion.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  
  if(results_common){
    
    d <- dmFit(d, dispersion = "common_dispersion", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
    
    d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
    
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
  d <- dmDispersion(d, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_mode = disp_mode, disp_init = common_disp, disp_moderation = disp_moderation, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  common_dispersion(d) <- common_disp
  
  plotDispersion(d, out_dir = out_name)
  
  d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  d <- dmTest(d, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  plotTest(d, out_dir = out_name)
  
  res <- results(d)
  
  write.table(res, paste0(out_name, "results.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  save(d, file = paste0(out_name, "d.Rdata"))
  
  rm("d")
  gc()
  
  
}


sessionInfo()







