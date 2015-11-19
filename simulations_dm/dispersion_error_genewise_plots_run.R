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

files <- list.files(out_dir, pattern = ".txt")
files



##############################################################################
### Plots
##############################################################################

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]

files_list <- paste0(out_dir, files)


for(i in 1:length(files)){
  # i = 1
  
  out_dir <- file_path_sans_ext(files_list[i])
  
  est <- read.table(files_list[i], header = TRUE, sep = "\t")
  
  est$dispersion <- factor(est$dispersion, levels = c("genewise", "moderated", "common"))
  est$method <- factor(est$method, levels = c("ML-dirmult", "PL", "CR"))
  
  
  ### Error as a difference
  
  error <- data.frame(error = est$est - est$true, dispersion = est$dispersion, method = est$method)
  error$dispersionmethod <- interaction(error$dispersion, error$method)
  
  ### Absolute error
  
  error$error <- abs(error$error)
  
  ylim <- c(min(aggregate(. ~ dispersionmethod, error[, c("error", "dispersionmethod")], whisker_lower)[, 2]) - 1, max(aggregate(. ~ dispersionmethod, error[, c("error", "dispersionmethod")], whisker_upper)[, 2]) + 1)
  
  
  ggp <- ggplot(data = error, aes(y = error, x = dispersion, fill = method)) + 
    theme_bw() +
    ylab("Absolute error") +
    geom_boxplot(outlier.size = 0) +
    coord_cartesian(ylim = ylim) +
    theme(axis.text = element_text(size = 16), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_boxplot_absolute.pdf"), 5, 5)
  print(ggp)
  dev.off()
  
  
  ggp <- ggplot(data = error, aes(y = log10(error), x = dispersion, fill = method)) + 
    theme_bw() +
    ylab("Log 10 of absolute error") +
    geom_boxplot(outlier.size = 1) +
    theme(axis.text = element_text(size = 16), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_boxplot_absolute_log.pdf"), 5, 5)
  print(ggp)
  dev.off()
  
  
  
  ### Error as a ratio
  
  error <- data.frame(error = est$est/est$true, dispersion = est$dispersion, method = est$method)
  error$dispersionmethod <- interaction(error$dispersion, error$method)
  
  
  ggp <- ggplot(data = error, aes(y = log10(error), x = dispersion, fill = method)) + 
    theme_bw() +
    ylab("Log10 of error ratio") +
    geom_boxplot(outlier.size = 1) +
    geom_hline(yintercept = 0, color="grey50", linetype = 2, size = 0.5) +
    theme(axis.text = element_text(size = 16), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_boxplot_ratio_log.pdf"), 5, 5)
  print(ggp)
  dev.off()
  
  #     ggp <- ggplot(data = error, aes(y = log10(error), x = dispersion, fill = method)) + 
  #       theme_bw() +
  #       ylab("Log10 of error ratio") +
  #       geom_violin(trim = FALSE, position = position_dodge(width = 0.9)) +
  #       geom_boxplot(outlier.size = 1, alpha = 0, position = position_dodge(width = 0.9), width = 0.1) +
  #       theme(axis.text = element_text(size = 16), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16))
  #     
  #     pdf(paste0(out_dir, "_violin_ratio_log.pdf"))
  #     print(ggp)
  #     dev.off()
  
  
}


















































