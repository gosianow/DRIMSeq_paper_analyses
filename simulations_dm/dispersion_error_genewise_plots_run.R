######################################################
## ----- dispersion_error_genewise_plots_run
## <<dispersion_error_genewise_plots_run.R>>

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

out_dir <- "error_genewise/"

error_files <- list.files(out_dir, pattern = ".txt")
error_files



##############################################################################
### Plots for common dispersion 
##############################################################################

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]

out_dir_list <- paste0(out_dir, error_files)


for(i in 1:length(error_files)){
  # i = 1
  
  out_dir <- file_path_sans_ext(out_dir_list[i])
  
  error <- read.table(out_dir_list[i], header = TRUE)
  
  error$dispersion <- factor(error$dispersion, levels = c("common", "genewise", "moderated"))
  error$method <- factor(error$method, levels = c("ML-dirmult", "PL", "CR"))
  
  colnames(error) <- c("Error", "Dispersion", "Method")
  error$DispersionMethod <- interaction(error$Dispersion, error$Method)
  
  
  ### Error
  
  ylim <- c(min(aggregate(. ~ DispersionMethod, error[, c("Error", "DispersionMethod")], whisker_lower)[, 2]) - 1, max(aggregate(. ~ DispersionMethod, error[, c("Error", "DispersionMethod")], whisker_upper)[, 2]) + 1)
  
  ggp <- ggplot(data = error, aes(y = Error, x = Dispersion, fill = Method)) + 
    theme_bw() +
    geom_boxplot(outlier.size = 0) +
    coord_cartesian(ylim = ylim) +
    geom_hline(yintercept = 0, color="grey", linetype = 2, size = 0.5) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  
  pdf(paste0(out_dir, "_boxplot_raw.pdf"))
  print(ggp)
  dev.off()
  
  
  ### Absolute error
  
  error$Error <- abs(error$Error)
  
  ylim <- c(min(aggregate(. ~ DispersionMethod, error[, c("Error", "DispersionMethod")], whisker_lower)[, 2]) - 1, max(aggregate(. ~ DispersionMethod, error[, c("Error", "DispersionMethod")], whisker_upper)[, 2]) + 1)
  
  
  ggp <- ggplot(data = error, aes(y = Error, x = Dispersion, fill = Method)) + 
    theme_bw() +
    ylab("Absolute error") +
    geom_boxplot(outlier.size = 0) +
    coord_cartesian(ylim = ylim) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_boxplot_absolute.pdf"))
  print(ggp)
  dev.off()
  
  
  ggp <- ggplot(data = error, aes(y = log10(Error), x = Dispersion, fill = Method)) + 
    theme_bw() +
    ylab("Log10 of absolute error") +
    geom_violin(trim = FALSE) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_violin_absolute_log.pdf"))
  print(ggp)
  dev.off()
  
  
  #   ggp <- ggplot(data = error, aes(y = Error, x = Dispersion, fill = Method)) + 
  #     theme_bw() +
  #     ylab("Absolute error") +
  #     geom_violin(trim = FALSE) +
  #     # coord_cartesian(ylim = ylim) +
  #     ylim(ylim[1], ylim[2]) +
  #     theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  # 
  #   pdf(paste0(out_dir, "_violin_absolute.pdf"))
  #   print(ggp)
  #   dev.off()
  
  
  ### Squared error
  
  error$Error <- error$Error^2
  
  ylim <- c(min(aggregate(. ~ DispersionMethod, error[, c("Error", "DispersionMethod")], whisker_lower)[, 2]) - 1, max(aggregate(. ~ DispersionMethod, error[, c("Error", "DispersionMethod")], whisker_upper)[, 2]) + 1)
  
  
  ggp <- ggplot(data = error, aes(y = Error, x = Dispersion, fill = Method)) + 
    theme_bw() +
    ylab("Squared error") +
    geom_boxplot(outlier.size = 0) +
    coord_cartesian(ylim = ylim) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  
  pdf(paste0(out_dir, "_boxplot_squared.pdf"))
  print(ggp)
  dev.off()
  
  
}


















































