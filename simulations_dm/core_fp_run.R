##############################################################################
## <<core_fp_run.R>>

# BioC 3.2
# Created 4 Mar 2016

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
# m=100 # Number of genes
# n=6 # Number of samples
# nm=100
# ### Common dispersion of gene expression
# nd=0
# disp_prior_df=0.1
# param_pi_path='/home/gosia/multinomial_project/simulations_dm/drimseq/dm_parameters_drimseq_0_3_3/kim_kallisto/prop_q3_kim_kallisto_fcutoff.txt'
# ### Genewise dispersion of feature proportions
# param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq/dm_parameters_drimseq_0_3_3/kim_kallisto/disp_common_kim_kallisto.txt'


##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(simulation_script)
print(workers)
print(sim_name)
print(run)
print(m)
print(n)
print(disp_prior_df)
print(nm)
print(nd)
print(param_pi_path)
print(param_gamma_path)



##############################################################################

source(simulation_script)


### Proportions
pi <- read.table(param_pi_path, header = FALSE, sep = "\t", as.is = TRUE)
pi <- as.numeric(pi[, 1])
pi <- pi/sum(pi) ### Make sure sum(pi) = 1

print(pi)

stopifnot(length(pi) > 1)

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
  g0_meanlog <- params[1, 2]
  g0_sdlog <- params[2, 2]
  print(g0_meanlog)
  print(g0_sdlog)
  sim_disp_genewise <- TRUE
}



##############################################################################

dir.create(rwd, recursive = T, showWarnings = FALSE)
setwd(rwd)

out_dir <- "core_fp/run/"
dir.create(out_dir, recursive = T, showWarnings = FALSE)

out_suffix <- "core_fp"


out_name <- paste0(sim_name, "n", n, "_nm", nm, "_nd", nd, "_", basename(file_path_sans_ext(param_pi_path)), "_",  basename(file_path_sans_ext(param_gamma_path)), "_")

out_name

if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}


##############################################################################
### Simulations for moderated dispersion 
##############################################################################

### Generate random dispersion
if(sim_disp_genewise)
  g0 <- rlnorm(m, meanlog = g0_meanlog, sdlog = g0_sdlog)


counts <- dm_simulate(m = m, n = 2*n, pi = pi, g0 = g0, nm = nm, nd = nd, mc.cores = workers)
print(head(counts))

group_split <- strsplit2(rownames(counts), ":")

d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep(c("c1", "c2"), each = n))


est <- list()
fp <- list()

##############################################################################
### No CR adjustement

d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = disp_prior_df, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)


d <- dmFit(d, dispersion = "common_dispersion", BPPARAM = BPPARAM)
d <- dmTest(d, BPPARAM = BPPARAM)
res <- results(d)

est[[1]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = common_dispersion(d), true = g0, pvalue = res$pvalue, dispersion = "common", method = "PL")
fp[[1]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), dispersion = "common", method = "PL")



d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
d <- dmTest(d, BPPARAM = BPPARAM)
res <- results(d)

est[[2]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = genewise_dispersion(d)$genewise_dispersion, true = g0, pvalue = res$pvalue, dispersion = "genewise", method = "PL")
fp[[2]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), dispersion = "genewise", method = "PL")



d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = common_dispersion(d), disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = disp_prior_df, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)

d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
d <- dmTest(d, BPPARAM = BPPARAM)
res <- results(d)

est[[3]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = genewise_dispersion(d)$genewise_dispersion, true = g0, pvalue = res$pvalue, dispersion = "moderated", method = "PL")
fp[[3]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), dispersion = "moderated", method = "PL")




##############################################################################
### With CR adjustement


d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = disp_prior_df, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)


d <- dmFit(d, dispersion = "common_dispersion", BPPARAM = BPPARAM)
d <- dmTest(d, BPPARAM = BPPARAM)
res <- results(d)

est[[4]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = common_dispersion(d), true = g0, pvalue = res$pvalue, dispersion = "common", method = "CR")
fp[[4]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), dispersion = "common", method = "CR")



d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
d <- dmTest(d, BPPARAM = BPPARAM)
res <- results(d)

est[[5]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = genewise_dispersion(d)$genewise_dispersion, true = g0, pvalue = res$pvalue, dispersion = "genewise", method = "CR")
fp[[5]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), dispersion = "genewise", method = "CR")



d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = common_dispersion(d), disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = disp_prior_df, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)


d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
d <- dmTest(d, BPPARAM = BPPARAM)
res <- results(d)

est[[6]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = genewise_dispersion(d)$genewise_dispersion, true = g0, pvalue = res$pvalue, dispersion = "moderated", method = "CR")
fp[[6]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), dispersion = "moderated", method = "CR")


##############################################################################
## True dispersion

d <- dmDispersion(d, mean_expression = TRUE, common_dispersion = FALSE, genewise_dispersion = FALSE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = disp_prior_df, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)

if(length(g0) == 1){
  gd <- rep(g0, m)
  names(gd) <- names(d@counts)
  genewise_dispersion(d) <- gd
}else{
  gd <- g0
  names(gd) <- names(d@counts)
  genewise_dispersion(d) <- gd
}


d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
d <- dmTest(d, BPPARAM = BPPARAM)
res <- results(d)

est[[7]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = genewise_dispersion(d)$genewise_dispersion, true = g0, pvalue = res$pvalue, dispersion = "genewise", method = "true")
fp[[7]] <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), dispersion = "genewise", method = "true")



##############################################################################
### Analysis with dirmult package

indx_group1 <- which(samples(d)$group == "c1")
indx_group2 <- which(samples(d)$group == "c2")


dirmult_results <- lapply(1:length(d@counts), function(g){
  # g = 1
  
  y0 <- t(d@counts[[g]])
  out0 <- dirmult(data = y0, trace = FALSE, epsilon = 10^(4))
  
  y1 <- t(d@counts[[g]][, indx_group1])
  out1 <- dirmult(data = y1, trace = FALSE, epsilon = 10^(4))
  
  y2 <- t(d@counts[[g]][, indx_group2])
  out2 <- dirmult(data = y2, trace = FALSE, epsilon = 10^(4))
  
  out_df <- data.frame(loglik0 = out0$loglik, loglik1 = out1$loglik, loglik2 = out2$loglik, gamma0 = sum(out0$gamma), gamma1 = sum(out1$gamma), gamma2 = sum(out2$gamma))
  
  return(out_df)
  
})

dirmult_results <- rbind.fill(dirmult_results)

lik_null <- dirmult_results$loglik0
lik_full <- dirmult_results$loglik1 + dirmult_results$loglik2

lr <- 2 * (lik_full - lik_null)
df <- results(d)$df
pvalue <- pchisq(lr, df = df , lower.tail = FALSE)


est[[8]] <- data.frame(gene_id = genewise_dispersion(d)$gene_id, est = dirmult_results$gamma0, true = g0, pvalue = pvalue, dispersion = "genewise", method = "ML-dirmult")
fp[[8]] <- data.frame(fp = mean(pvalue < 0.05, na.rm = TRUE), dispersion = "genewise", method = "ML-dirmult")


est <- rbind.fill(est)
fp <- rbind.fill(fp)



write.table(est, paste0(out_dir, out_name, "est_", out_suffix, "_", run,".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(fp, paste0(out_dir, out_name, "fp_", out_suffix, "_", run,".txt"), quote = FALSE, sep = "\t", row.names = FALSE)


sessionInfo()


