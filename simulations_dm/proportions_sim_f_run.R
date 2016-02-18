######################################################
## ----- proportions_sim_f_run
## <<proportions_sim_f_run.R>>

# BioC 3.2 - R32dev
# Created 10 Feb 2016 


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
# m=1000
# n=3 # Number of samples
# nm=1000
# nd=0
# nr_features=c(3:10)
# param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq/dm_parameters_drimseq_0_3_3/kim_kallisto/disp_common_kim_kallisto.txt'
# save=TRUE
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



##############################################################################

dir.create(rwd, recursive = T, showWarnings = FALSE)
setwd(rwd)

out_dir <- "proportions_f/run/"
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
  # j = 3
  print(nr_features[j])
  
  switch(out_suffix, 
         
         proportions_decay = {
           ### Proportions with decay
           pi <- (1/2)^(1:nr_features[j])
           
         },
         
         proportions_uniform = {
           ### Proportions uniform
           pi <- rep(1, nr_features[j])
           
         }
  )
  
  pi <- sort(pi/sum(pi), decreasing = TRUE) ### Make sure sum(pi) = 1
  print(pi)
  
  counts <- dm_simulate(m = m, n = 2*n, pi = pi, g0 = g0, nm = nm, nd = nd, mc.cores = workers)
  print(head(counts))
  
  group_split <- strsplit2(rownames(counts), ":")
  
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep(c("c1", "c2"), each = n))
  
  
  d <- dmFilter(d, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0, max_features = Inf)
  
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)
  
  d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
  
  ## LR test
  d <- dmTest(d, test = "lr", BPPARAM = BPPARAM)
  
  res <- results(d)
  
  est[[paste0("moderation_none_lr", j)]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = round(genewise_dispersion(d)$genewise_dispersion, 2), true = round(g0, 2), pvalue = res$pvalue, test_statistic = res[, "lr"], df = res$df, nr_features = nr_features[j], disp_estimator = "moderation_none", test = "lr")
  
  fp[[paste0("moderation_none_lr", j)]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), nr_features = nr_features[j], disp_estimator = "moderation_none", test = "lr")
  
  
  ## F test
  d <- dmTest(d, test = "f", BPPARAM = BPPARAM)
  
  res <- results(d)
  
  est[[paste0("moderation_none_f", j)]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = round(genewise_dispersion(d)$genewise_dispersion, 2), true = round(g0, 2), pvalue = res$pvalue, test_statistic = res[, "f"], df1 = res$df1, df2 = res$df2, nr_features = nr_features[j], disp_estimator = "moderation_none", test = "f")
  
  fp[[paste0("moderation_none_f", j)]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), nr_features = nr_features[j], disp_estimator = "moderation_none", test = "f")
  
  
  if(save)
    save("d", file = paste0(out_dir, out_name, "d_", out_suffix, "_", run, "_", nr_features[j], "_moderation_none",".Rdata"))
  
  ### True dispersion 
  
  common_dispersion(d) <- g0
  
  d <- dmFit(d, dispersion = "common_dispersion", BPPARAM = BPPARAM)
  
  ## LR test
  d <- dmTest(d, test = "lr", BPPARAM = BPPARAM)
  
  res <- results(d)
  
  est[[paste0("true_lr", j)]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = round(g0, 2), true = round(g0, 2), pvalue = res$pvalue, test_statistic = res[, "lr"], df = res$df, nr_features = nr_features[j], disp_estimator = "true", test = "lr")
  
  fp[[paste0("true_lr", j)]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), nr_features = nr_features[j], disp_estimator = "true", test = "lr")
  
  
  ## F test
  d <- dmTest(d, test = "f", BPPARAM = BPPARAM)
  
  res <- results(d)
  
  est[[paste0("true_f", j)]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = round(g0, 2), true = round(g0, 2), pvalue = res$pvalue, test_statistic = res[, "f"], df1 = res$df1, df2 = res$df2, nr_features = nr_features[j], disp_estimator = "true", test = "f")
  
  fp[[paste0("true_f", j)]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), nr_features = nr_features[j], disp_estimator = "true", test = "f")
  
  
  if(save)
    save("d", file = paste0(out_dir, out_name, "d_", out_suffix, "_", run, "_", nr_features[j], "_true",".Rdata"))
  
  rm("d")
  
}


est <- rbind.fill(est)
fp <- rbind.fill(fp)


write.table(est, paste0(out_dir, out_name, "est_", out_suffix, "_", run,".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(fp, paste0(out_dir, out_name, "fp_", out_suffix, "_", run,".txt"), quote = FALSE, sep = "\t", row.names = FALSE)


sessionInfo()









