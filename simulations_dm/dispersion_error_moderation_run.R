######################################################
## ----- dispersion_error_moderation_run
## <<dispersion_error_moderation_run.R>>

# BioC 3.1
# Created 6 Nov 2015 

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

# rwd='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/'
# workers=2
# sim_name=''
# r=1 # Number of simulations
# m=100 # Number of genes
# n=3 # Number of samples
# nm=100 # Mean gene expression
# nd=0 # Negative binomial dispersion of gene expression
# disp_prior_df=seq(0,1,by=1)
# param_pi_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/prop_q3_uniform.txt'
# param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/disp_genewise_kim_kallisto_lognormal.txt'


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

out_dir <- "error_moderation/"
dir.create(out_dir, recursive = T, showWarnings = FALSE)

out_dir <- paste0("error_moderation/", sim_name, "n", n, "_nm", nm, "_nd", nd, "_", basename(file_path_sans_ext(param_pi_path)), "_",  basename(file_path_sans_ext(param_gamma_path)), "_")


##############################################################################
### Simulations for moderated dispersion 
##############################################################################


error_list <- lapply(1:r, function(i){
  # i = 1
  
  # g0 <- rgamma(m, shape = g0_shape, scale = g0_scale)
  g0 <- rlnorm(m, meanlog = g0_meanlog, sdlog = g0_sdlog)
  gamma_est <- data.frame(true = g0)
  
  counts <- dm_simulate(m = m, n = n, pi = pi, g0 = g0, nm = nm, nd = nd, mc.cores = workers)
  print(head(counts))
  
  group_split <- strsplit2(rownames(counts), ":")
  
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep("c1", ncol(counts)))
  
  disp_interval <- c(0, max(g0)*100)
  
  ### With CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = FALSE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  common_disp <- common_dispersion(d)
  
  ### Dispersion with different degree of moderation 
  
  mse <- list()
  for(j in 1:length(disp_prior_df)){
    # j = 1
    
    d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = common_disp, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = disp_prior_df[j], disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
    common_dispersion(d) <- common_disp
    
    gamma_est[, paste0(disp_prior_df[j])] <- genewise_dispersion(d)$genewise_dispersion
    
    mse[[j]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, disp_prior_df = disp_prior_df[j])
    
  }
  
  ### Plot the distribution of estimated gamma_+
  
  gamma_estm <- melt(gamma_est, variable.name = "Moderation", value.name = "Dispersion")
  
  pdf(paste0(out_dir, "gammas.pdf"))
  ggp <- ggplot(gamma_estm, aes(x = Moderation, y = log10(Dispersion))) + 
    geom_violin(trim = FALSE) + 
    ylab("Log 10 of dispersion gamma_+") +
    geom_hline(yintercept = log10(common_disp), color = "blue")
  print(ggp)
  dev.off()
  
  mse <- do.call(rbind, mse)
  
  print(mem_used())
  rm("d")
  gc()
  print(mem_used())
  
  return(mse)
  
})




error <- do.call(rbind, error_list)

write.table(error, paste0(out_dir, "error_moderation.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


##############################################################################
### Plots for moderated dispersion 
##############################################################################














sessionInfo()


