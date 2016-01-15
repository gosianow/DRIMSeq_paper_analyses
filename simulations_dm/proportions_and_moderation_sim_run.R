######################################################
## ----- proportions_and_moderation_sim_run
## <<proportions_and_moderation_sim_run.R>>

# BioC 3.2
# Created 14 Jan 2015 


##############################################################################
Sys.time()
##############################################################################

library(BiocParallel)
library(pryr)
library(plyr)
library(dirmult)
library(limma)
library(DRIMSeq)
library(ggplot2)
library(reshape2)
library(tools)


##############################################################################
# Arguments for testing the code
##############################################################################

# rwd='/home/gosia/multinomial_project/simulations_dm/drimseq/'
# simulation_script='/home/gosia/R/drimseq_paper/simulations_dm/dm_simulate.R'
# workers=4
# sim_name='test_'
# run='run1'
# m=100
# n=3 # Number of samples
# nm=1000
# nd=0
# nr_features=c(2,4)
# param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq/dm_parameters_drimseq_0_3_3/kim_kallisto/disp_genewise_kim_kallisto_lognormal.txt'
# out_suffix='proportions_decay'


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
print(simulation_script)
print(workers)
print(sim_name)
print(m)
print(n)
print(nm)
print(nd)
print(nr_features)
print(param_gamma_path)

##############################################################################


source(simulation_script)



### Dispersion
params <- read.table(param_gamma_path, header = FALSE, sep = "\t")

sim_disp_common <- FALSE
sim_disp_genewise <- FALSE

# Common dispersion
if(ncol(params) == 1){
  g0 <- as.numeric(params)
  print(g0)
  sim_disp_common <- TRUE
}

# Genewise dispersion from lognormal distribution
if(ncol(params) == 2){
  g0_meanlog <- params[1, 2]
  g0_sdlog <- params[2, 2]
  print(g0_meanlog)
  print(g0_sdlog)
  sim_disp_genewise <- TRUE
}



stopifnot(any(c("proportions_decay", "proportions_uniform") %in% out_suffix))

pi_function <- function(nr_features, out_suffix){
  
  if(out_suffix == "proportions_decay"){
    
    pi <- (1/2)^(1:nr_features)
    
  }else if(out_suffix == "proportions_uniform"){
    
    pi <- rep(1, nr_features)
    
  }
  
  pi <- sort(pi/sum(pi), decreasing = TRUE) ### Make sure sum(pi) = 1
  pi
}



##############################################################################

dir.create(rwd, recursive = T, showWarnings = FALSE)
setwd(rwd)

out_dir <- "proportions_and_moderation/run/"
dir.create(out_dir, recursive = T, showWarnings = FALSE)


out_name <- paste0(sim_name, "n", n, "_nm", nm, "_nd", nd, "_",  basename(file_path_sans_ext(param_gamma_path)), "_")

out_name


if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}



##############################################################################
### Simulations from proportions with decay 
##############################################################################


est <- list()
fp <- list()


for(j in 1:length(nr_features)){
  # j = 1
  print(nr_features[j])
  
  ### Random dispersion
  if(sim_disp_genewise)
    g0 <- rlnorm(m, meanlog = g0_meanlog, sdlog = g0_sdlog)

  ### Proportions
  pi <- pi_function(nr_features[j], out_suffix)
  print(pi)
  
  counts <- dm_simulate(m = m, n = 2*n, pi = pi, g0 = g0, nm = nm, nd = nd, mc.cores = workers)
  print(head(counts))
  
  group_split <- strsplit2(rownames(counts), ":")
  
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep(c("c1", "c2"), each = n))
  
  
  d_tmp <- dmFilter(d, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0, max_features = Inf)
  
  
  ### Use no moderation for dispersion estimation
  
  d <- dmDispersion(d_tmp, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 0.1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)
  
  common_disp <- common_dispersion(d)
  
  d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
  d <- dmTest(d, BPPARAM = BPPARAM)
  res <- results(d)
  
  est[[paste0("moderation_none", j)]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = round(genewise_dispersion(d)$genewise_dispersion, 2), true = round(g0, 2), pvalue = res$pvalue, nr_features = nr_features[j], disp_estimator = "moderation_none")
  
  fp[[paste0("moderation_none", j)]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), nr_features = nr_features[j], disp_estimator = "moderation_none")
  
  rm("d")
  
  ### Use common moderation for dispersion estimation
  
  d <- dmDispersion(d_tmp, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = common_disp, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = 0.1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)
  
  
  d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
  d <- dmTest(d, BPPARAM = BPPARAM)
  res <- results(d)
  
  est[[paste0("moderation_common", j)]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = round(genewise_dispersion(d)$genewise_dispersion, 2), true = round(g0, 2), pvalue = res$pvalue, nr_features = nr_features[j], disp_estimator = "moderation_common")
  
  fp[[paste0("moderation_common", j)]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), nr_features = nr_features[j], disp_estimator = "moderation_common")
  
  
  ### Use true dispersion estimates
  names(g0) <- genewise_dispersion(d)$gene_id
  genewise_dispersion(d) <- g0
  
  d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
  d <- dmTest(d, BPPARAM = BPPARAM)
  res <- results(d)
  
  est[[paste0("true", j)]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = round(g0, 2), true = round(g0, 2), pvalue = res$pvalue, nr_features = nr_features[j], disp_estimator = "true")
  
  fp[[paste0("true", j)]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), nr_features = nr_features[j], disp_estimator = "true")
  
  rm("d")
  rm("d_tmp")
  
}


est <- rbind.fill(est)
fp <- rbind.fill(fp)


write.table(est, paste0(out_dir, out_name, "est_", out_suffix, "_", run,".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(fp, paste0(out_dir, out_name, "fp_", out_suffix, "_", run,".txt"), quote = FALSE, sep = "\t", row.names = FALSE)



sessionInfo()




