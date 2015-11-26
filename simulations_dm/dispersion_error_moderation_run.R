######################################################
## ----- dispersion_error_moderation_run
## <<dispersion_error_moderation_run.R>>

# BioC 3.1
# Created 6 Nov 2015 

##############################################################################
library(BiocParallel)
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
workers=4
sim_name=''
r=1 # Number of simulations
m=100 # Number of genes
n=3 # Number of samples
nm=100 # Mean gene expression
nd=0 # Negative binomial dispersion of gene expression
disp_prior_df=seq(0,1,by=1)
param_pi_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/prop_q3_uniform.txt'
### Genewise dispersion
# param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/disp_genewise_kim_kallisto_lognormal.txt'
### Common dispersion
param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/disp_common_kim_kallisto.txt'

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

### Proportions
pi <- as.numeric(read.table(param_pi_path, header = FALSE, sep = "\t", as.is = TRUE)[, 1])
pi <- pi/sum(pi) ### Make sure sum(pi) = 1
print(pi)

### Dispersion
params <- read.table(param_gamma_path, header = FALSE, sep = "\t")

sim_disp_common <- FALSE
sim_disp_genewise <- FALSE

### Common dispersion
if(ncol(params) == 1){
  g0 <- as.numeric(params)
  print(g0)
  sim_disp_common <- TRUE
}

### Genewise dispersion from lognormal distribution
if(ncol(params) == 2){
  params <- as.numeric(params[, 2])
  g0_meanlog <- params[1]
  g0_sdlog <- params[2]
  print(g0_meanlog)
  print(g0_sdlog)
  sim_disp_genewise <- TRUE
}



dir.create(rwd, recursive = T, showWarnings = FALSE)
setwd(rwd)

out_dir <- "error_moderation/"
dir.create(out_dir, recursive = T, showWarnings = FALSE)

out_dir <- paste0("error_moderation/", sim_name, "n", n, "_nm", nm, "_nd", nd, "_", basename(file_path_sans_ext(param_pi_path)), "_",  basename(file_path_sans_ext(param_gamma_path)), "_")


if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
  }else{
    BPPARAM <- SerialParam()
  }


##############################################################################
### Simulations for moderated dispersion 
##############################################################################


est_list <- lapply(1:r, function(i){
  # i = 1
  print(i)
  
  if(sim_disp_genewise)
    g0 <- rlnorm(m, meanlog = g0_meanlog, sdlog = g0_sdlog)
  
  
  counts <- dm_simulate(m = m, n = n, pi = pi, g0 = g0, nm = nm, nd = nd, mc.cores = workers)
  print(head(counts))
  
  group_split <- strsplit2(rownames(counts), ":")
  
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep("c1", ncol(counts)))
  
  ### With CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = FALSE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)
  
  common_disp <- common_dispersion(d)
  
  ### Dispersion with different degree of moderation 
  
  est <- list()
  
  for(j in 1:length(disp_prior_df)){
    # j = 1
    
    d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = common_disp, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = disp_prior_df[j], disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)
    common_dispersion(d) <- common_disp
    
    est[[j]] <- data.frame(est = genewise_dispersion(d)$genewise_dispersion, true = g0, disp_prior_df = disp_prior_df[j])
    
  }
  
  est <- do.call(rbind, est)
  
  rm("d")
  gc()
  
  return(est)
  
})




est <- do.call(rbind, est_list)

write.table(est, paste0(out_dir, "est_moderation.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


sessionInfo()


