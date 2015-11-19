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

files <- list.files(out_dir, pattern = ".txt")
files


##############################################################################
### Plots for moderated dispersion 
##############################################################################

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]

files_list <- paste0(out_dir, files)


for(i in 1:length(files)){
  # i = 1
  
  out_dir <- file_path_sans_ext(files_list[i])
  
  est <- read.table(files_list[i], header = TRUE)
  
  est$disp_prior_df <- factor(est$disp_prior_df)
  
  ### Error as a difference
  
  error <- data.frame(error = est$est - est$true, disp_prior_df = est$disp_prior_df)
  
  ### Absolute error
  
  error$error <- abs(error$error)
  
  min_median <- min(aggregate(. ~ disp_prior_df, error, median)[, 2])
  ylim <- c(min(aggregate(. ~ disp_prior_df, error, whisker_lower)[, 2]) - 1, max(aggregate(. ~ disp_prior_df, error, whisker_upper)[, 2]) + 1)
  

#   ggp <- ggplot(data = error, aes(y = error, x = disp_prior_df)) + 
#     theme_bw() +
#     ylab("Absolute error") +
#     xlab("Moderation") +
#     geom_violin(trim = FALSE, show_guide = FALSE, fill = "grey80", colour = "grey80") +
#     geom_boxplot(outlier.size = 1, fill = NA, width = 0.1, alpha = 0.1) +
#     coord_cartesian(ylim = ylim) +
#     geom_hline(yintercept = min_median, color = "red", linetype = 2, size = 0.3) +
#     theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
#   
#   pdf(paste0(out_dir, "_violin_absolute.pdf"))
#   print(ggp)
#   dev.off()
  
  
  
  ggp <- ggplot(data = error, aes(y = log10(error), x = disp_prior_df)) + 
    theme_bw() +
    ylab("Log10 of absolute error") +
    xlab("Moderation") +
    geom_violin(trim = FALSE, show_guide = FALSE, fill = "grey80", colour = "grey80") +
    geom_boxplot(outlier.size = 1.2, fill = NA, width = 0.5) +
    geom_hline(yintercept = log10(min_median), color="red", linetype = 2, size = 0.5) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_violin_absolute_log.pdf"))
  print(ggp)
  dev.off()
  
  
  
  ### Error as a ratio
  
  error <- data.frame(error = est$est/est$true, disp_prior_df = est$disp_prior_df)

  ggp <- ggplot(data = error, aes(y = log10(error), x = disp_prior_df)) + 
    theme_bw() +
    ylab("Log10 of error ratio") +
    xlab("Moderation") +
    geom_violin(trim = FALSE, show_guide = FALSE, fill = "grey80", colour = "grey80") +
    geom_boxplot(outlier.size = 1.2, fill = NA, width = 0.5) +
    geom_hline(yintercept = 0, color="red", linetype = 2, size = 0.5) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_violin_ratio_log.pdf"))
  print(ggp)
  dev.off()
  
  
  
  
}



sessionInfo()


