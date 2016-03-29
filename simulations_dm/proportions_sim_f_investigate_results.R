
# R32dev

# library(devtools)
# load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")

# source("/home/gosia/R/drimseq_paper/simulations_dm/proportions_sim_f_investigate_results.R")

##############################################################################
# Simulations from DM - proportions_f
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

rwd='/home/gosia/multinomial_project/simulations_dm/drimseq/'
simulation_script='/home/gosia/R/drimseq_paper/simulations_dm/dm_simulate.R'
workers=5
sim_name=''
run='run1'
m=1000
n=3 # Number of samples
nm=1000
nd=0
nr_features=c(3,8)
param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq/dm_parameters_drimseq_0_3_3/kim_kallisto/disp_common_kim_kallisto.txt'
save=FALSE
out_suffix='proportions_uniform'

##############################################################################
# Read in the arguments
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


plot_dir <- paste0("proportions_f/investigate_results/", out_suffix, "/")
dir.create(plot_dir, recursive = T, showWarnings = FALSE)


##############################################################################
### Load simulation results 
##############################################################################



for(j in 1:length(nr_features)){
  # j = 1
  
  print(nr_features[j])
  
  out_name <- paste0(sim_name, "n", n, "_nm", nm, "_nd", nd, "_",  basename(file_path_sans_ext(param_gamma_path)), "_")
  
  load(paste0(out_dir, out_name, "d_", out_suffix, "_", run, "_", nr_features[j], "_true",".Rdata"))
  
  dt <- d  
  
  load(paste0(out_dir, out_name, "d_", out_suffix, "_", run, "_", nr_features[j], "_moderation_none",".Rdata"))
  
  dm <- d
  
  out_name <- paste0(plot_dir, out_name, nr_features[j], "_")
  
  
  ## LR test ; true dispersion
  dt <- dmTest(dt, test = "lr", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  rest <- results(dt)
  plotTest(dt, out_dir = paste0(out_name, "dt_"))
  
  
  ## LR test ; estimated - moderated dispersion
  dm <- dmTest(dm, test = "lr", BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  resm <- results(dm)
  plotTest(dm, out_dir = paste0(out_name, "dm_"))
  
  
  
  ### Plot genes with p-value ~ 0 for estimated dispersion (and corresponding true)
  # 
  # table(resm$adj_pvalue < 0.05)
  # 
  # oo <- order(resm$pvalue, decreasing = FALSE)
  # genes <- resm[oo, "gene_id"]
  # 
  # 
  # for(i in 1:10){
  #   plotFit(dm, gene_id = genes[i], out_dir = paste0(out_name, "dm", "_xdm0_", i, "_"))
  # }
  # 
  # 
  # for(i in 1:10){
  #   plotFit(dt, gene_id = genes[i], out_dir = paste0(out_name, "dt", "_xdm0_", i, "_"))
  # }
  
  
  
  
  
  ### Plot dipsersion versus p-values
  
  resmm <- merge(resm, genewise_dispersion(dm), by = "gene_id", sort = FALSE)
  
  resmm <- resmm[order(resmm$pvalue, decreasing = FALSE), ]
  
  
  ggp <- ggplot(data = resmm, aes(x = pvalue, y = genewise_dispersion)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = common_dispersion(dt), slope = 0, color = "orange") +
    stat_smooth()
  
  
  pdf(paste0(out_name, "dm", "_disp_vs_pvalues.pdf"))
  print(ggp)
  dev.off()
  
  
  
  
  
  ### Plot difference between group dispersions versus p-values
  
  dm1 <- dmDispersion(dm[ ,samples(dm)$group == "c1"], mean_expression = FALSE,
    common_dispersion = FALSE, genewise_dispersion = TRUE,
    disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
    disp_tol = 1e-08, disp_init = common_dispersion(dm), disp_init_weirMoM = TRUE,
    disp_grid_length = 21, disp_grid_range = c(-10, 10),
    disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
    prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0,
    BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  
  dm2 <- dmDispersion(dm[ ,samples(dm)$group == "c2"], mean_expression = FALSE,
    common_dispersion = FALSE, genewise_dispersion = TRUE,
    disp_adjust = TRUE, disp_mode = "grid", disp_interval = c(0, 1e+05),
    disp_tol = 1e-08, disp_init = common_dispersion(dm), disp_init_weirMoM = TRUE,
    disp_grid_length = 21, disp_grid_range = c(-10, 10),
    disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3,
    prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = 0,
    BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  
  resm12m <- merge(genewise_dispersion(dm1), genewise_dispersion(dm2), by = "gene_id", sort = FALSE, suffixes = c("1", "2"))
  
  resm12m <- merge(resm, resm12m, by = "gene_id", sort = FALSE)
  resm12m$disp_diff <- abs(resm12m$genewise_dispersion1 - resm12m$genewise_dispersion2)
  
  ggp <- ggplot(data = resm12m, aes(x = pvalue, y = log10(disp_diff))) +
    geom_point(alpha = 0.6) +
    stat_smooth()
  
  
  
  pdf(paste0(out_name, "dm", "_disp_diff_vs_pvalues.pdf"))
  print(ggp)
  dev.off()
  
  
  
  
  
  
  
  
  ### Plot p-values versus p-values with true and estimated dispersion 
  
  true_disp <- common_dispersion(dt)
  
  resmtm <- merge(resm, rest, by = "gene_id", sort = FALSE, suffixes = c("m", "t"))
  resmtm <- merge(resmtm, genewise_dispersion(dm), by = "gene_id", sort = FALSE)
  
  ggp <- ggplot(data = resmtm, aes(x = log10(pvaluet), y = log10(pvaluem), colour = genewise_dispersion)) +
    # geom_bin2d(binwidth = c(0.05, 0.05)) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "black") +
    geom_hline(yintercept = log10(0.05), color = "blue") +
    geom_vline(xintercept = log10(0.05), color = "blue") +
    scale_colour_gradient2(low = "brown", mid = "orange", high = "darkgreen", midpoint = true_disp)
  
  
  pdf(paste0(out_name, "dm", "_pval_vs_pval.pdf"), 9, 7)
  print(ggp)
  dev.off()
  
  
  table(resmtm$genewise_dispersion < true_disp & resmtm$pvaluem < resmtm$pvaluet)
  
  
  
  
  ### Plot the difference in proportion estimates
  
  prop_diff <- proportions(dt)[, -c(1,2)] - proportions(dm)[, -c(1,2)]
  
  prop_diff <- data.frame(proportions(dt)[, c(1,2)], prop_diff)
  
  prop_diff <- merge(prop_diff, genewise_dispersion(dm), by = "gene_id", sort = FALSE)
  
  prop_diff$feature_id <- factor(prop_diff$feature_id, levels = unique(prop_diff$feature_id))
  
  
  ggp <- ggplot(data = prop_diff, aes(x = genewise_dispersion, y = c1)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, color = "blue") +
    facet_wrap(~ feature_id, nrow = 2, scales = "free_y")
  
  pdf(paste0(out_name, "dm", "_prop_diff_c1.pdf"), 14, 7)
  print(ggp)
  dev.off()
  
  ggp <- ggplot(data = prop_diff, aes(x = genewise_dispersion, y = c2)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, color = "blue") +
    facet_wrap(~ feature_id, nrow = 2, scales = "free_y")
  
  pdf(paste0(out_name, "dm", "_prop_diff_c2.pdf"), 14, 7)
  print(ggp)
  dev.off()
  
  ggp <- ggplot(data = prop_diff, aes(x = genewise_dispersion, y = null)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, color = "blue") +
    facet_wrap(~ feature_id, nrow = 2, scales = "free_y")
  
  pdf(paste0(out_name, "dm", "_prop_diff_null.pdf"), 14, 7)
  print(ggp)
  dev.off()
  
  
  
  ### Plot genes with pvalue < 0.05 for estimated dispersion and pvalue > 0.05 for true dispersion
  
  resmtm <- merge(resm, rest, by = "gene_id", sort = FALSE, suffixes = c("m", "t"))
  
  resmtm <- resmtm[order(resmtm$pvaluem, decreasing = FALSE), ]
  resmtm <- resmtm[order(resmtm$pvaluet, decreasing = TRUE), ]
  
  resmtm <- resmtm[resmtm$pvaluem < 0.05, ]
  
  genes <- resmtm$gene_id
  
  
  
  for(i in 1:10){
    plotFit(dm, gene_id = genes[i], out_dir = paste0(out_name, "dm", "_xdm0dt1_", i, "_"))
  }
  
  for(i in 1:10){
    plotFit(dt, gene_id = genes[i], out_dir = paste0(out_name, "dt", "_xdm0dt1_", i, "_"))
  }
  
  
  
  
  
  ### Check the contribution of each transcript to the LR test
  
  lik_contribution <- function(y, pi, gamma0){
    
    l <- rowSums(lgamma(y + pi * gamma0) - lgamma(pi * gamma0), na.rm = TRUE)
    
    return(l)
    
  }
  
  
  
  # d = dt; name = "dt"
  
  lr_contribution <- function(d, name){
    
    lrc <- lapply(1:length(d), function(i){
      # i = 2
      gamma0 <- slot(d, d@dispersion)
      
      if(length(gamma0) > 1)
        gamma0 <- gamma0[i]
      
      ## Null
      y <- d@counts[[i]]
      pi <- d@fit_null[[i]][, "null"]
      
      likn <- lik_contribution(y, pi, gamma0)
      
      ### Check if calculation is correct
      N <- ncol(y)
      S <- colSums(y)
      l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + sum(likn, na.rm = TRUE)
      
      d@fit_null@metadata[i, "lik"] - l
      
      
      ### Full
      groups <- levels(d@samples$group)
      
      likf <- sapply(1:length(groups), function(ii){
        # ii = 1
        
        group <- groups[ii]
        y <- d@counts[[i]][, d@samples$group == group]
        pi <- d@fit_full[[i]][, group]
        
        likfg <- lik_contribution(y, pi, gamma0)
        
        ### Check if calculation is correct
        N <- ncol(y)
        S <- colSums(y)
        l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + sum(likfg, na.rm = TRUE)
        
        d@fit_full@metadata[i, paste0("lik_", group)] - l
        
        return(likfg)
        
      })
      
      lrc <- 2 * (rowSums(likf, na.rm = TRUE) - likn)
      
      data.frame(t(lrc))
      
    })
    
    
    lrc <- rbind.fill(lrc)
    
    lrc$gene_id <- names(d)
    
    lrcm <- melt(lrc, id.vars = "gene_id", variable.name = "feature", value.name = "lr_contrib")
    
    lrcm <- merge(lrcm, results(d), by = "gene_id", sort = FALSE)
    
    lrcm$fp <- lrcm$pvalue < 0.05
    
    
    ggp <- ggplot(data = lrcm, aes(x = feature, y = lr_contrib, colour = feature)) + 
      geom_jitter(alpha = 0.3, show.legend = FALSE) +
      geom_boxplot(colour = "black", fill = NA, outlier.size = NA) +
      theme_bw()
    
    
    pdf(paste0(out_name, name, "_lr_contribution.pdf"), 10, 7)
    print(ggp)
    dev.off()
    
    
    
    ggp <- ggplot(data = lrcm, aes(x = feature, y = lr_contrib, colour = fp, fill = fp)) + 
      geom_jitter(alpha = 0.3, show.legend = FALSE, position = position_jitterdodge(dodge.width = 0.75)) +
      geom_boxplot(fill = NA, outlier.size = NA) +
      theme_bw()
    
    
    pdf(paste0(out_name, name, "_lr_contribution_fp.pdf"), 10, 7)
    print(ggp)
    dev.off()
    
    return(invisible(lrcm))
    
  }
  
  lrc_dt <- lr_contribution(d = dt, name = "dt")
  
  lrc_dm <- lr_contribution(d = dm, name = "dm")
  
  
  
  
  ### Plot LR contribution for genes with pvalue < 0.05 for estimated dispersion and pvalue > 0.05 for true dispersion
  
  resmtm <- merge(resm, rest, by = "gene_id", sort = FALSE, suffixes = c("m", "t"))
  
  resmtm <- resmtm[order(resmtm$pvaluem, decreasing = FALSE), ]
  resmtm <- resmtm[order(resmtm$pvaluet, decreasing = TRUE), ]
  
  resmtm <- resmtm[resmtm$pvaluem < 0.05, ]
  
  genes <- resmtm$gene_id
  
  
  for(i in 1:10){
    # i = 1
    
    lrc_tmp <- lrc_dm[lrc_dm$gene_id == genes[i], , drop = FALSE]
    
    lr_tmp <- sum(lrc_tmp$lr_contrib, na.rm = TRUE)
    
    print(all.equal(lr_tmp, lrc_tmp$lr[1]))
    
    ggp <- ggplot(lrc_tmp, aes(x = feature, y = lr_contrib)) +
      geom_bar(stat = "identity") +
      ggtitle(paste0("LR = ", lr_tmp))
    
    
    pdf(paste0(out_name, "dm", "_xdm0dt1_", i, "_", "lr_contribution_",genes[i],".pdf"))
    print(ggp)
    dev.off()
    
  }
  
  for(i in 1:10){
    # i = 1
    
    lrc_tmp <- lrc_dt[lrc_dt$gene_id == genes[i], , drop = FALSE]
    
    lr_tmp <- sum(lrc_tmp$lr_contrib, na.rm = TRUE)
    
    print(all.equal(lr_tmp, lrc_tmp$lr[1]))
    
    ggp <- ggplot(lrc_tmp, aes(x = feature, y = lr_contrib)) +
      geom_bar(stat = "identity") +
      ggtitle(paste0("LR = ", lr_tmp))
    
    
    pdf(paste0(out_name, "dt", "_xdm0dt1_", i, "_", "lr_contribution_",genes[i],".pdf"))
    print(ggp)
    dev.off()
    
  }
  
  
  
  
}


























