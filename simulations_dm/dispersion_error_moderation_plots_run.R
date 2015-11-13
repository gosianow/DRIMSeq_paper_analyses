######################################################
## ----- dispersion_error_moderation_plots_run
## <<dispersion_error_moderation_plots_run.R>>

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

out_dir <- "error_moderation/"

error_files <- list.files(out_dir, pattern = ".txt")



##############################################################################
### Plots for moderated dispersion 
##############################################################################

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]

error_files_list <- paste0(out_dir, error_files)


for(i in 1:length(error_files)){
  # i = 1
  
  error <- read.table(error_files_list[i], header = TRUE)
  
  out_dir <- file_path_sans_ext(error_files_list[i])
  
  error$disp_prior_df <- factor(error$disp_prior_df)
  colnames(error) <- c("Error", "Moderation")
  
  ### Error
  
  min_median <- min(aggregate(. ~ Moderation, error, median)[, 2])
  
  ylim <- c(min(aggregate(. ~ Moderation, error, whisker_lower)[, 2]) - 1, max(aggregate(. ~ Moderation, error, whisker_upper)[, 2]) + 1)
  
  ggp <- ggplot(data = error, aes(y = Error, x = Moderation)) + 
    theme_bw() +
    ylab("Error") +
    geom_boxplot(outlier.size = 0, fill = "grey80") +
    coord_cartesian(ylim = ylim) +
    geom_hline(yintercept = 0, color = "grey", linetype = 1, size = 0.5) +
    geom_hline(yintercept = min_median, color = "red", linetype = 2, size = 0.3) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  
  pdf(paste0(out_dir, "_boxplot_raw.pdf"))
  print(ggp)
  dev.off()
  
  
  ### Absolute error
  
  error$Error <- abs(error$Error)
  
  min_median <- min(aggregate(. ~ Moderation, error, median)[, 2])
  
  ylim <- c(min(aggregate(. ~ Moderation, error, whisker_lower)[, 2]) - 1, max(aggregate(. ~ Moderation, error, whisker_upper)[, 2]) + 1)
  
  ggp <- ggplot(data = error, aes(y = Error, x = Moderation)) + 
    theme_bw() +
    ylab("Absolute error") +
    geom_boxplot(outlier.size = 0, fill = "grey80") +
    coord_cartesian(ylim = ylim) +
    geom_hline(yintercept = min_median, color = "red", linetype = 2, size = 0.3) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_boxplot_absolute.pdf"))
  print(ggp)
  dev.off()
  
  
  ggp <- ggplot(data = error, aes(y = Error, x = Moderation)) + 
    theme_bw() +
    ylab("Absolute error") +
    geom_violin(trim = FALSE, show_guide = FALSE, fill = "grey80", colour = "grey80") +
    geom_boxplot(outlier.size = 1, fill = NA, width = 0.1, alpha = 0.1) +
    coord_cartesian(ylim = ylim) +
    geom_hline(yintercept = min_median, color = "red", linetype = 2, size = 0.3) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_violin_absolute.pdf"))
  print(ggp)
  dev.off()
  
  
  
  # log
  ggp <- ggplot(data = error, aes(y = log10(Error), x = Moderation)) + 
    theme_bw() +
    ylab("Log 10 of absolute error") +
    geom_boxplot(outlier.size = 1, fill = "grey80") +
    geom_hline(yintercept = log10(min_median), color="red", linetype = 2, size = 0.3) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  
  pdf(paste0(out_dir, "_boxplot_absolute_log.pdf"))
  print(ggp)
  dev.off()
  
  
  
  ggp <- ggplot(data = error, aes(y = log10(Error), x = Moderation)) + 
    theme_bw() +
    ylab("Log10 of absolute error") +
    geom_violin(trim = FALSE, show_guide = FALSE, fill = "grey80", colour = "grey80") +
    geom_boxplot(outlier.size = 1.2, fill = "grey80", width = 0.5) +
    geom_hline(yintercept = log10(min_median), color="red", linetype = 2, size = 0.3) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_violin_absolute_log.pdf"))
  print(ggp)
  dev.off()
  
  
  
  ### Squared error
  
  error$Error <- error$Error^2
  
  min_median <- min(aggregate(. ~ Moderation, error, median)[, 2])
  
  ylim <- c(min(aggregate(. ~ Moderation, error, whisker_lower)[, 2]) - 1, max(aggregate(. ~ Moderation, error, whisker_upper)[, 2]) + 1)
  
  ggp <- ggplot(data = error, aes(y = Error, x = Moderation)) + 
    theme_bw() +
    ylab("Squared error") +
    geom_boxplot(outlier.size = 0, fill = "grey80") +
    coord_cartesian(ylim = ylim) +
    geom_hline(yintercept = min_median, color = "red", linetype = 2, size = 0.3) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  
  pdf(paste0(out_dir, "_boxplot_squared.pdf"))
  print(ggp)
  dev.off()
  
  
  # log
  ggp <- ggplot(data = error, aes(y = log10(Error), x = Moderation)) + 
    theme_bw() +
    ylab("Log 10 of squared error") +
    geom_boxplot(outlier.size = 1, fill = "grey80") +
    geom_hline(yintercept = log10(min_median), color="red", linetype = 2, size = 0.3) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  
  pdf(paste0(out_dir, "_boxplot_squared_log.pdf"))
  print(ggp)
  dev.off()

  
  ### violinplot
  
  # ggp <- ggplot(data = error, aes(y = Error, x = Moderation)) + 
  #   theme_bw() +
  #   ylab("Squared error") +
  #   geom_violin(trim = FALSE, fill = "grey80") +
  #   theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  # 
  # pdf(paste0(out_dir, "_violin_squared.pdf"))
  # print(ggp)
  # dev.off()
  
  # log
  
  # ggp <- ggplot(data = error, aes(y = log10(Error), x = Moderation)) + 
  #   theme_bw() +
  #   ylab("Log 10 of squared error") +
  #   geom_violin(trim = FALSE, fill = "grey80") +
  #   theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  # 
  # pdf(paste0(out_dir, "_violin_squared_log.pdf"))
  # print(ggp)
  # dev.off()
  
  
  
  
}



sessionInfo()


