######################################################
## ----- moderation_real_auto_run
## <<moderation_real_auto_run.R>>

# BioC 3.2
# Created 16 Apr 2016 

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

rwd='/home/gosia/multinomial_project/simulations_dm/drimseq'
out_dir='moderation_real_auto_v2/run'
simulation_script='/home/gosia/R/drimseq_paper/simulations_dm/dm_simulate.R'
workers=5
sim_name=''
run='run1'
m=200
n=3
param_pi_path='/home/gosia/multinomial_project/simulations_dm/drimseq/dm_parameters_drimseq_0_3_3/brooks_kallisto/prop_brooks_kallisto_fcutoff.txt'
param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq/dm_parameters_drimseq_0_3_3/brooks_kallisto/disp_genewise_brooks_kallisto_lognormal.txt'
param_nm_path='/home/gosia/multinomial_project/simulations_dm/drimseq/dm_parameters_drimseq_0_3_3/brooks_kallisto/nm_brooks_kallisto_lognormal.txt'
param_nd_path='/home/gosia/multinomial_project/simulations_dm/drimseq/dm_parameters_drimseq_0_3_3/brooks_kallisto/nd_common_brooks_kallisto.txt'

##############################################################################
# Read in the arguments
##############################################################################

rm(list = ls())

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

cat(paste0(args, collapse = "\n"), fill = TRUE)


print(rwd)
print(simulation_script)
print(workers)
print(sim_name)
print(run)
print(m)
print(n)
print(param_nm_path)
print(param_nd_path)
print(param_pi_path)
print(param_gamma_path)


##############################################################################

source(simulation_script)


### Proportions
pi <- read.table(param_pi_path, header = TRUE, sep = "\t", as.is = TRUE)
pi_list <- split(pi$proportions, pi$gene_id)


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




### Mean gene expression
params <- read.table(param_nm_path, header = FALSE, sep = "\t")

nm_meanlog <- params[1, 2]
nm_sdlog <- params[2, 2]

print(nm_meanlog)
print(nm_sdlog)


# Negative binomial dispersion of gene expression
params <- read.table(param_nd_path, header = FALSE, sep = "\t")

nd <- as.numeric(params)

print(nd)


##############################################################################

dir.create(rwd, recursive = T, showWarnings = FALSE)
setwd(rwd)


dir.create(out_dir, recursive = T, showWarnings = FALSE)

out_suffix <- "moderation_real_auto"


out_name <- paste0(sim_name, "n", n, "_", basename(file_path_sans_ext(param_nm_path)), "_", basename(file_path_sans_ext(param_nd_path)), "_", basename(file_path_sans_ext(param_pi_path)), "_",  basename(file_path_sans_ext(param_gamma_path)), "_")

out_name

if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}


##############################################################################
### Simulations for moderated dispersion 
##############################################################################



### Random dispersion
if(sim_disp_genewise)
  g0 <- rlnorm(m, meanlog = g0_meanlog, sdlog = g0_sdlog)

### Random proportions
pi <- pi_list[sample(1:length(pi_list), m, replace = TRUE)]

### Random gene expression
nm <- round(rlnorm(m, meanlog = nm_meanlog, sdlog = nm_sdlog))

counts <- dm_simulate(m = m, n = 2*n, pi = pi, g0 = g0, nm = nm, nd = nd, mc.cores = workers)
print(head(counts))

group_split <- strsplit2(rownames(counts), ":")

d_org <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep(c("c1", "c2"), each = n))

names(g0) <- names(d_org@counts)
names(pi) <- names(d_org@counts)
names(nm) <- names(d_org@counts)

d <- dmFilter(d_org, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0, max_features = Inf)

keep_genes <- names(d@counts)




### With CR adjustement

d_org <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = FALSE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.1, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)

common_disp <- common_dispersion(d_org)



### Dispersion with automated degree of moderation 


d <- dmDispersion(d_org, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05), disp_tol = 1e-08, disp_init = common_disp, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = 0.1, disp_span = 0.1, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BPPARAM)

disp_prior_df <- d@priorn

common_dispersion(d) <- common_disp


est <- data.frame(est = round(genewise_dispersion(d)$genewise_dispersion, 2), true = round(g0[keep_genes], 2), disp_prior_df = disp_prior_df, q = sapply(pi[keep_genes], length), nm = nm[keep_genes])


d <- dmFit(d, dispersion = "genewise_dispersion", BPPARAM = BPPARAM)
d <- dmTest(d, BPPARAM = BPPARAM)
res <- results(d)

fp <- data.frame(fp = mean(res$pvalue < 0.05, na.rm = TRUE), disp_prior_df = disp_prior_df)
est$pvalue <- res$pvalue






write.table(est, paste0(out_dir, "/", out_name, "est_", out_suffix, "_", run,".txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(fp, paste0(out_dir, "/", out_name, "fp_", out_suffix, "_", run,".txt"), quote = FALSE, sep = "\t", row.names = FALSE)







sessionInfo()


