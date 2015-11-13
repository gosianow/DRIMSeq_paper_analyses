######################################################
## ----- dispersion_error_common_run
## <<dispersion_error_common_run.R>>

# BioC 3.1
# Created 11 Nov 2015 

##############################################################################

library(pryr)
library(dirmult)
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
workers=2
sim_name=''
r=1 # Number of simulations
m=100 # Number of genes
n=3 # Number of samples
nm=100 # Mean gene expression
nd=0 # Negative binomial dispersion of gene expression
disp_prior_df=0.1 # One value
param_pi_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/prop_q3_uniform.txt'
param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/disp_common_kim_kallisto.txt' # One value


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
print(param_gamma_path)
print(disp_prior_df)


##############################################################################

pi <- as.numeric(read.table(param_pi_path, header = FALSE, sep = "\t", as.is = TRUE)[, 1])
pi <- pi/sum(pi) ### Make sure sum(pi) = 1

g0 <- as.numeric(read.table(param_gamma_path, header = FALSE, sep = "\t"))

print(pi)
print(g0)

# # Proportions
# pi <- c(1/3, 1/3, 1/3) 
# # Common dispersion gamma_+
# g0 <- 28

dir.create(rwd, recursive = T, showWarnings = FALSE)
setwd(rwd)

out_dir <- "error_common/"
dir.create(out_dir, recursive = T, showWarnings = FALSE)

out_dir <- paste0("error_common/", sim_name, "n", n, "_nm", nm, "_nd", nd, "_", basename(file_path_sans_ext(param_pi_path)), "_",  basename(file_path_sans_ext(param_gamma_path)), "_")


##############################################################################
### Simulations for common dispersion 
##############################################################################


error_list <- lapply(1:r, function(i){
  # i = 1
  
  mse <- list()
  
  counts <- dm_simulate(m = m, n = n, pi = pi, g0 = g0, nm = nm, nd = nd, mc.cores = workers)
  print(head(counts))
  
  group_split <- strsplit2(rownames(counts), ":")
  
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep("c1", ncol(counts)))
  
  disp_interval <- c(0, max(g0)*100)
  
  ### No CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[1]] <- data.frame(mse = common_dispersion(d) - g0, dispersion = "common", method = "PL")
  
  mse[[2]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "genewise", method = "PL")
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = common_dispersion(d), disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[3]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "moderated", method = "PL")
  
  
  ### With CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[4]] <- data.frame(mse = common_dispersion(d) - g0, dispersion = "common", method = "CR")
  
  mse[[5]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "genewise", method = "CR")
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = common_dispersion(d), disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[6]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "moderated", method = "CR")
  
  
  ### gammas estimated with dirmult package
  
  dirmult_dispersion <- unlist(lapply(1:length(d@counts), function(g){
    # g = 1
    y <- t(d@counts[[g]])
    out <- dirmult(data = y, trace = FALSE, epsilon = 10^(4))
    sum(out$gamma)
    # (1-out$theta)/out$theta
  }))
  
  mse[[7]] <- data.frame(mse = dirmult_dispersion - g0, dispersion = "genewise", method = "ML-dirmult")
  
  mse <- do.call(rbind, mse)
  return(mse)
  
})




error <- do.call(rbind, error_list)

write.table(error, paste0(out_dir, "error_common.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


##############################################################################
### Plots for common dispersion 
##############################################################################





sessionInfo()












