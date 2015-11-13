######################################################
## ----- dispersion_error_optimization_plots_run
## <<dispersion_error_optimization_plots_run.R>>

# BioC 3.1
# Created 12 Nov 2015 

##############################################################################

library(ggplot2)
library(reshape2)
library(tools)


##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/'


##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)


##############################################################################

setwd(rwd)

out_dir <- "error_optimization/"

error_files <- list.files(out_dir, pattern = "error_optimization.txt")
error_files

estimates_files <- list.files(out_dir, pattern = "estimates_optimization.txt")
estimates_files


##############################################################################
### Plots for genewise dispersion - different optimization methods 
##############################################################################

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]


error_files_list <- paste0(out_dir, error_files)
estimates_files_list <- paste0(out_dir, estimates_files)


for(i in 1:length(error_files)){
  # i = 1
  
  ###### Plot errors 
  
  out_dir <- file_path_sans_ext(error_files_list[i])
  
  error <- read.table(error_files_list[i], header = TRUE)
  
  error$method <- factor(error$method, levels = c("grid", "optimize", "optim", "constrOptim"))
  colnames(error) <- c("Error", "Method")
  
  
  ### Error
  
  ylim <- c(min(aggregate(. ~ Method, error, whisker_lower)[, 2]) - 1, max(aggregate(. ~ Method, error, whisker_upper)[, 2]) + 1)
  
  ggp <- ggplot(data = error, aes(y = Error, x = Method, fill = Method)) + 
    theme_bw() +
    geom_boxplot(outlier.size = 0, show_guide = FALSE) +
    coord_cartesian(ylim = ylim) +
    geom_hline(yintercept = 0, color="grey", linetype = 2, size = 0.5) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  
  pdf(paste0(out_dir, "_boxplot_raw.pdf"))
  print(ggp)
  dev.off()
  
  
    ggp <- ggplot(data = error, aes(y = Error, x = Method, fill = Method)) + 
      theme_bw() +
      ylab("Error") +
      geom_violin(trim = FALSE, show_guide = FALSE) +
      # coord_cartesian(ylim = ylim) +
      # ylim(ylim[1], ylim[2]) +
      theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
    pdf(paste0(out_dir, "_violin_raw.pdf"))
    print(ggp)
    dev.off()
  
  
  ### Absolute error
  
  error$Error <- abs(error$Error)
  
  ylim <- c(min(aggregate(. ~ Method, error, whisker_lower)[, 2]) - 1, max(aggregate(. ~ Method, error, whisker_upper)[, 2]) + 1)
  
  ggp <- ggplot(data = error, aes(y = Error, x = Method, fill = Method)) + 
    theme_bw() +
    ylab("Absolute error") +
    geom_boxplot(outlier.size = 0, show_guide = FALSE) +
    coord_cartesian(ylim = ylim) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_boxplot_absolute.pdf"))
  print(ggp)
  dev.off()
  
  
  ggp <- ggplot(data = error, aes(y = log10(Error), x = Method, fill = Method)) + 
    theme_bw() +
    ylab("Log 10 of absolute error") +
    geom_boxplot(outlier.size = 1, show_guide = FALSE) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_boxplot_absolute_log.pdf"))
  print(ggp)
  dev.off()
  
  
  ggp <- ggplot(data = error, aes(y = log10(Error), x = Method, fill = Method)) + 
    theme_bw() +
    ylab("Log10 of absolute error") +
    geom_violin(trim = FALSE, show_guide = FALSE) +
    # coord_cartesian(ylim = ylim) +
    # ylim(ylim[1], ylim[2]) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_violin_absolute_log.pdf"))
  print(ggp)
  dev.off()
  
  
  ### Squared error
  
  error$Error <- error$Error^2
  
  ylim <- c(min(aggregate(. ~ Method, error, whisker_lower)[, 2]) - 1, max(aggregate(. ~ Method, error, whisker_upper)[, 2]) + 1)
  
  ggp <- ggplot(data = error, aes(y = Error, x = Method, fill = Method)) + 
    theme_bw() +
    ylab("Squared error") +
    geom_boxplot(outlier.size = 0, show_guide = FALSE) +
    coord_cartesian(ylim = ylim) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  
  pdf(paste0(out_dir, "_boxplot_squared.pdf"))
  print(ggp)
  dev.off()
  
  
  ###### Plot estimates

  out_dir <- file_path_sans_ext(estimates_files_list[i])
  
  estimates <- read.table(estimates_files_list[i], header = TRUE)
  
  true_disp <- estimates[1, "true"]
  
  estimatesm <- melt(estimates[, c("grid", "optimize", "optim", "constrOptim")])
  estimatesm$variable <- factor(estimatesm$variable, levels = c("grid", "optimize", "optim", "constrOptim"))
  levels(estimatesm$variable)
  
  
  ggp <- ggplot(data = estimatesm, aes(y = log10(value), x = variable, fill = variable)) + 
    theme_bw() +
    xlab("Method") +
    ylab("Log 10 of gamma_+") +
    geom_violin(trim = FALSE, show_guide = FALSE) +
    geom_hline(yintercept = log10(true_disp), color="black", linetype = 2, size = 0.5) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_violin_log.pdf"))
  print(ggp)
  dev.off()
  
  
}














































