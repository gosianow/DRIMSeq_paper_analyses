######################################################
## ----- dispersion_error_diff_gamma_run
## <<dispersion_error_diff_gamma_run.R>>

# BioC 3.1
# Created 28 Nov 2015 

##############################################################################

library(BiocParallel)
library(pryr)
library(limma)
library(DRIMSeq)
library(ggplot2)
library(reshape2)
library(tools)

source("/home/gosia/R/drimseq_paper/simulations_dm/dm_simulate.R")

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/'
workers=4
sim_name=''
r=1 # Number of simulations
m=500 # Number of genes
n=3 # Number of samples
nm=100 # Mean gene expression
nd=0 # Negative binomial dispersion of gene expression
param_pi_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/prop_q10_kim_kallisto_overall.txt'
g0=c(1,5,10,20,30,100,200,1000)

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
print(r)
print(m)
print(n)
print(nm)
print(nd)
print(sim_name)
print(param_pi_path)
print(g0)


##############################################################################

pi <- as.numeric(read.table(param_pi_path, header = FALSE, sep = "\t", as.is = TRUE)[, 1])
pi <- pi/sum(pi) ### Make sure sum(pi) = 1

print(pi)
print(g0)


dir.create(rwd, recursive = T, showWarnings = FALSE)
setwd(rwd)

out_dir <- "error_diff_gamma/"
dir.create(out_dir, recursive = T, showWarnings = FALSE)

out_name <- paste0(out_dir, sim_name, "n", n, "_nm", nm, "_nd", nd, "_", basename(file_path_sans_ext(param_pi_path)), "_")


if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}


##############################################################################
### Simulations for genewise dispersion - different optimization methods 
##############################################################################



est_list <- lapply(1:r, function(i){
  # i = 1
  print(i)
  
  est <- list()
  
  for(j in 1:length(g0)){
    
    print(j)
    
    counts <- dm_simulate(m = m, n = n, pi = pi, g0 = g0[j], nm = nm, nd = nd, mc.cores = workers)
    print(head(counts))
    
    group_split <- strsplit2(rownames(counts), ":")
    
    d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep("c1", ncol(counts)))
    
    
    ## grid
    
    d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)
    
    est[[j]] <- data.frame(est = genewise_dispersion(d)$genewise_dispersion, true = g0[j])
    
    
  }
  
  
  est <- do.call(rbind, est)
  rownames(est) <- NULL
  
  rm("d")
  gc()
  
  return(est)
  
})



est <- do.call(rbind, est_list)

write.table(est, paste0(out_name, "est_diff_gamma.txt"), quote = FALSE, sep = "\t", row.names = FALSE)



sessionInfo()






