######################################################
## ----- sim5_drimseq_0_3_1_run
## <<sim5_drimseq_0_3_1_run.R>>

# BioC 3.1
# Created 16 Nov 2015 

##############################################################################

library(DRIMSeq)
library(limma)
library(ggplot2)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_sim5'
simulation='drosophila_node_nonull'
workers=4
count_method=c('htseq','kallisto','htseq_prefiltered15')[2]
dispersion_common=TRUE
results_common=TRUE
disp_mode_list=c('grid','grid','optimize','optim','constrOptim')[1]
disp_moderation_list=c('none','common','none','none','none')[1]


##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(simulation)
print(workers)
print(count_method)
print(dispersion_common)
print(results_common)
print(disp_mode_list)
print(disp_moderation_list)


##############################################################################

setwd(paste0(rwd, "/", simulation))
method_out <- "drimseq_0_3_1"

########################################################
# create metadata file
########################################################

metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata


##########################################################################
# Load counts
##########################################################################

out_dir <- paste0(method_out, "/", count_method, "/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


if(count_method == "htseq")
  count_dir <- "2_counts/dexseq_nomerge/dexseq"
if(count_method == "kallisto")
  count_dir <- "2_counts/kallisto/kallisto"
if(count_method == "htseq_prefiltered15")
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast15/dexseq"


### load counts
counts_list <- lapply(1:6, function(i){
  # i = 1
  cts <- read.table(paste0(count_dir, i, ".txt"), header = FALSE, as.is = TRUE)
  colnames(cts) <- c("group_id", paste0("sample_", i))  
  return(cts)
})

counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)

counts <- counts[!grepl(pattern = "_", counts$group_id),]
group_split <- strsplit2(counts[,1], ":")
counts <- counts[, -1]



##########################################################################
# DRIMSeq analysis
##########################################################################


d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sample_id, group = metadata$group)


### Filtering
table(samples(d)$group)
d <- dmFilter(d, min_samps_gene_expr = 3, min_samps_feature_prop = 3, min_feature_prop = 0.01)


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







