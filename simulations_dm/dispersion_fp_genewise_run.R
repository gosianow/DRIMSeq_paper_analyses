######################################################
## ----- dispersion_fp_genewise_run
## <<dispersion_fp_genewise_run.R>>

# BioC 3.1
# Created 16 Nov 2015 

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
workers=5
sim_name='test_'
r=3 # Number of simulations
m=100 # Number of genes
n=2 # Number of samples
nm=100 # Mean gene expression
nd=0 # Negative binomial dispersion of gene expression
disp_prior_df=0.1 # One value
param_pi_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/prop_q3_uniform.txt'
param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/disp_genewise_kim_kallisto_lognormal.txt' # Parameters for a distribution


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

params <- as.numeric(read.table(param_gamma_path, header = FALSE, sep = "\t")[, 2])

g0_meanlog <- params[1]
g0_sdlog <- params[2]

print(pi)
print(g0_meanlog)
print(g0_sdlog)

# # Proportions
# pi <- c(1/3, 1/3, 1/3) 
# # Genewise dispersion gamma_+
# g0_meanlog <- 3
# g0_sdlog <- 1



dir.create(rwd, recursive = T, showWarnings = FALSE)
setwd(rwd)

out_dir <- "fp_genewise/"
dir.create(out_dir, recursive = T, showWarnings = FALSE)

out_dir <- paste0("fp_genewise/", sim_name, "n", n, "_nm", nm, "_nd", nd, "_", basename(file_path_sans_ext(param_pi_path)), "_",  basename(file_path_sans_ext(param_gamma_path)), "_")

if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
  }else{
    BPPARAM <- SerialParam()
  }
  
##############################################################################
### Simulations for genewise dispersion 
##############################################################################


fp_list <- lapply(1:r, function(i){
  # i = 1
  print(i)
  
  fp <- list()
  
  g0 <- rlnorm(m, meanlog = g0_meanlog, sdlog = g0_sdlog)
  
  counts <- dm_simulate(m = m, n = n * 2, pi = pi, g0 = g0, nm = nm, nd = nd, mc.cores = workers)
  print(head(counts))
  
  group_split <- strsplit2(rownames(counts), ":")
  
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep(c("c1", "c2"), each = n))
  

  ### With CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)
  d <- dmFit(d, dispersion = "common_dispersion", BPPARAM = BPPARAM)
  d <- dmTest(d, BPPARAM = BPPARAM)
  res <- results(d)
  
  fp[[1]] <- data.frame(fp = mean(res$pvalue < 0.05), dispersion = "common", method = "CR")
  
  d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
  d <- dmTest(d, BPPARAM = BPPARAM)
  res <- results(d)
  
  fp[[2]] <- data.frame(fp = mean(res$pvalue < 0.05), dispersion = "genewise", method = "CR")
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = common_dispersion(d), disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)
  
  d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
  d <- dmTest(d, BPPARAM = BPPARAM)
  res <- results(d)
  
  fp[[3]] <- data.frame(fp = mean(res$pvalue < 0.05), dispersion = "moderated", method = "CR")
  
  fp <- do.call(rbind, fp)
  rownames(fp) <- NULL
  
  rm("d")
  gc()

  return(fp)
  
})




fp <- do.call(rbind, fp_list)
colnames(fp) <- c("fp", "dispersion", "method")

write.table(fp, paste0(out_dir, "fp_genewise.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


##############################################################################
### Plots for genewise dispersion 
##############################################################################





sessionInfo()
























