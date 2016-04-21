######################################################
## <<kim_drimseq_0_3_3_auto_moderation.R>>

# BioC 3.2
# Created 18 Apr 2016

##############################################################################
Sys.time()
##############################################################################

library(BiocParallel)
library(DRIMSeq)
library(limma)

##############################################################################
# Test arguments
##############################################################################

# rwd="/home/Shared/data/seq/kim_adenocarcinoma/"
# workers=5
# count_method=c("htseq", "kallisto")[1]
# model=c("model_full", "model_null_normal1", "model_null_tumor1")[2]
# method_out='drimseq_0_3_3'
# dmDS_auto_moderation_diagnostics_function_path='/home/gosia/R/drimseq_paper/help_functions/dmDS_auto_moderation_diagnostics.R'

##############################################################################
# Read in the arguments
##############################################################################


## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(workers)
print(count_method)
print(model)



##############################################################################

setwd(rwd)

if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}


out_dir <- paste0(method_out, "/",  model, "/", count_method, "/")
dir.create(out_dir, recursive = TRUE)

out_dir_tmp <- paste0(out_dir, "auto_moderation/")
dir.create(out_dir_tmp, recursive = TRUE)


##########################################################################
# Load DRIMSeq data
##########################################################################

### Load object d
load(paste0(out_dir, "drimseq_genewise_grid_none_d.Rdata"))



common_disp <- as.numeric(read.table(paste0(out_dir, "common_dispersion.txt")))
common_disp



###########################################################################
### Automatic moderation
###########################################################################

source(dmDS_auto_moderation_diagnostics_function_path)

x <- d

dmDS_auto_moderation_diagnostics(x = x, common_disp = common_disp, out_dir_tmp = out_dir_tmp, BPPARAM = BPPARAM)
























