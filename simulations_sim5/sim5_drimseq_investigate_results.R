##############################################################################
# Simulation sim5
##############################################################################

# R32

library(DRIMSeq)
library(limma)
library(ggplot2)
library(reshape2)

library(devtools)
load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")

##############################################################################
# Test arguments
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_sim5'
simulation='hsapiens_node_nonull'
workers=10
count_method=c('htseq','kallisto','htseqprefiltered15','htseqprefiltered5','kallistofiltered5','kallistoprefiltered5')[2]
filter_method=c("filter0", "filter2")[1]
dispersion_common=TRUE
results_common=TRUE
disp_mode=c('grid','grid','optimize','optim','constrOptim')[1]
disp_moderation=c('none','common','none','none','none')[1]
lr_contribution_function_path <- "/home/gosia/R/drimseq_paper/help_functions/dmDS_lr_contribution.R"

##############################################################################

setwd(paste0(rwd, "/", simulation))
method_out <- "drimseq_0_3_3"


source(lr_contribution_function_path)

########################################################
# create metadata file
########################################################

metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata


##########################################################################
# 
##########################################################################

out_dir <- paste0(method_out, "/", count_method, "/", filter_method, "/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



##########################################################################
# DRIMSeq results
##########################################################################


### Load object d

common_disp <- as.numeric(read.table(paste0(out_dir, "common_dispersion.txt")))
common_disp

disp <- "genewise"

if(disp_mode == "grid"){
  out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_", disp_moderation, "_")
}else{
  out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_")
}


load(paste0(out_name, "d.Rdata"))



### New out directory 
out_dir <- paste0(method_out, "/", count_method, "/", filter_method, "/", "investigate_results/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if(disp_mode == "grid"){
  out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_", disp_moderation, "_")
}else{
  out_name <- paste0(out_dir, "/drimseq_", disp, "_", disp_mode, "_")
}





res <- results(d)

plotTest(d, out_dir = paste0(out_name))




#######################################################
# load simulation info
#######################################################


simulation_details <- read.table("3_truth/simulation_details.txt", header = TRUE, as.is = TRUE, sep = "\t")

truth_file <- list.files("3_truth/", pattern = "truth")
truth_file

truth <- read.table(paste0("3_truth/", truth_file), header = TRUE, as.is = TRUE, sep = "\t")

truth <- truth[, c("gene", "ds_status", "de_status", "TPM", "nbr_isoforms", "diff_IsoPct", "nbrexonbins")]
rownames(truth) <- truth$gene

colnames(truth)[1] <- "gene_id"


##########################################################################
### transcript expected counts versus observed 

## observed counts
counts_ob <- counts(d)

## expected counts
counts_ex <- simulation_details[, c("gene_id", "transcript_id", colnames(simulation_details)[grepl("_isoformCount", colnames(simulation_details))])]

colnames(counts_ex)[c(1,2)] <- c("gene_id", "feature_id") 
colnames(counts_ex)[c(3:8)] <- paste0(colnames(counts_ob)[c(3:8)])

counts_m <- merge(counts_ob, counts_ex, by = c("gene_id", "feature_id"), sort = FALSE, suffixes = c("_ob","_ex"))


##########################################################################
### merge results with truth


resm <- merge(res, truth, by = c("gene_id"), sort = FALSE)

resm <- unique(resm)




##########################################################################
### plot proportions for FP ####


res_fp <- resm[resm$ds_status == 0 & resm$adj_pvalue < 0.05, ]

res_fp <- unique(res_fp)

res_fp <- res_fp[order(res_fp$pvalue, decreasing = FALSE), ]

genes <- res_fp[, "gene_id"]

for(i in 1:10){
  plotFit(d, gene_id = genes[i], out_dir = paste0(out_name, "xfp_", i, "_"), order = FALSE)
}



### Use expected counts
for(i in 1:10){
  # i= 1 
  counts_tmp <- counts_m[counts_m$gene_id == genes[i], ]
  
  dex <- dmDSdata(counts = counts_tmp[, grepl("_ex", colnames(counts_tmp))], gene_id = counts_tmp$gene_id, feature_id = counts_tmp$feature_id, sample_id = paste0("s", 1:6), group = rep(c("c1", "c2"), each = 3))
  
  
  dex <- dmDispersion(dex, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_init = common_dispersion(d))
  
  dex <- dmFit(dex)
  
  dex <- dmTest(dex)
  
  plotFit(dex, gene_id = genes[i], out_dir = paste0(out_name, "xfp_counts_ex_", i, "_"), order = FALSE)
  
}



### Plot the LR contribution for FP


for(i in 1:10){
  
  lr_contribution(d, gene_id = genes[i], name = paste0(out_name, "xfp_", i, "_"))
  
}



##########################################################################
### plot proportions for TP ----

res_tp <- resm[resm$ds_status == 1 & resm$adj_pvalue < 0.05, ]

res_tp <- res_tp[order(res_tp$pvalue, decreasing = FALSE), ]

genes <- res_tp[, "gene_id"]

for(i in 1:10){
  plotFit(d, gene_id = genes[i], out_dir = paste0(out_name, "xtp_", i, "_"), order = FALSE)
}



### Plot the LR contribution for TP

for(i in 1:10){
  # i = 2
  lr_contribution(d, gene_id = genes[i], name = paste0(out_name, "xtp_", i, "_"))
  
}





### Save an example data for Mark 

# save(d, gene_id, resm, file = "/home/gosia/tmp/example_lr_contributions.Rdata")
# 
# ### all results
# 
# head(resm)
# 
# ### plot LR contributions
# lr_contribution(d, gene_id, name = "")
# 
# ### plot proportions
# plotFit(d, gene_id, out_dir = "", order = FALSE)





##########################################################################
### Plot dispersion versus mean with marked FP ----

common_disp <- common_dispersion(d)

ggp <- plotDispersion(d)

df <- ggp$data

res_fp <- resm[resm$ds_status == 0 & resm$adj_pvalue < 0.05, ]

df_fp <- df[res_fp[, "gene_id"], ]

ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion")) +
  theme_bw() +
  xlab("Log10 of mean expression") +
  ylab("Log10 of gamma_+") +
  geom_point(size = 1.2, alpha = 0.4, na.rm = TRUE) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") + 
  geom_hline(yintercept = log10(common_disp), colour = "black", linetype = "dashed", size =  0.5) +
  geom_point(data = df_fp, aes_string(x = "mean_expression", y = "dispersion"), color = "orange")


pdf(paste0(out_name, "dispersion_vs_mean_fp.pdf"))
print(ggp)
dev.off()


##########################################################################
### Plot dispersion versus mean with marked TP ####

common_disp <- common_dispersion(d)

ggp <- plotDispersion(d)

df <- ggp$data

res_tp <- resm[resm$ds_status == 1 & resm$adj_pvalue < 0.05, ]

df_tp <- df[res_tp[, "gene_id"], ]

ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion")) +
  theme_bw() +
  xlab("Log10 of mean expression") +
  ylab("Log10 of gamma_+") +
  geom_point(size = 1.2, alpha = 0.4, na.rm = TRUE) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") + 
  geom_hline(yintercept = log10(common_disp), colour = "black", linetype = "dashed", size =  0.5) +
  geom_point(data = df_tp, aes_string(x = "mean_expression", y = "dispersion"), color = "blue")


pdf(paste0(out_name, "dispersion_vs_mean_tp.pdf"))
print(ggp)
dev.off()

##########################################################################
### Plot dispersion versus mean with marked P ####

common_disp <- common_dispersion(d)

ggp <- plotDispersion(d)

df <- ggp$data

res_p <- resm[resm$ds_status == 1, ]

df_p <- df[res_p[, "gene_id"], ]

ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion")) +
  theme_bw() +
  xlab("Log10 of mean expression") +
  ylab("Log10 of gamma_+") +
  geom_point(size = 1.2, alpha = 0.4, na.rm = TRUE) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") + 
  geom_hline(yintercept = log10(common_disp), colour = "black", linetype = "dashed", size =  0.5) +
  geom_point(data = df_p, aes_string(x = "mean_expression", y = "dispersion"), color = "dodgerblue")


pdf(paste0(out_name, "dispersion_vs_mean_p.pdf"))
print(ggp)
dev.off()





























##########################################################################
### Estimate dispersion per group ####


disp_init <- common_dispersion(dfql)
disp_init

dsub1 <- dfql[gene, 1:3]

dsub1 <- dmDispersion(dsub1, mean_expression = TRUE,
  common_dispersion = FALSE, genewise_dispersion = TRUE,
  disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
  disp_tol = 1e-12, disp_init = disp_init, disp_init_weirMoM = TRUE,
  disp_grid_length = 21, disp_grid_range = c(-10, 10),
  disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

genewise_dispersion(dsub1)


dsub1 <- dmFit(dsub1, dispersion = "genewise_dispersion",
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotFit(dsub1, gene_id = gene, out_dir = paste0(out_name, "dsub1", "_neg_dev_no_zeros_", genenr ,"_"))



dsub2 <- dfql[gene, 4:6]
samples(dsub2)


dsub2 <- dmDispersion(dsub2, mean_expression = TRUE,
  common_dispersion = FALSE, genewise_dispersion = TRUE,
  disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
  disp_tol = 1e-12, disp_init = disp_init, disp_init_weirMoM = TRUE,
  disp_grid_length = 21, disp_grid_range = c(-10, 10),
  disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

genewise_dispersion(dsub2)


dsub2 <- dmFit(dsub2, dispersion = "genewise_dispersion",
  prop_mode = "constrOptimG", prop_tol = 1e-16, verbose = 0,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers))

plotFit(dsub2, gene_id = gene, out_dir = paste0(out_name, "dsub2", "_neg_dev_no_zeros_", genenr ,"_"))

















































