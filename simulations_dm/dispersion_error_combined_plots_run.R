######################################################
## ----- dispersion_error_combined_plots_run
## <<dispersion_error_combined_plots_run.R>>

# BioC 3.1
# Created 17 Nov 2015 

##############################################################################

library(ggplot2)
library(reshape2)
library(tools)
library(limma)

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

out_dir_main <- "error_combined/"
dir.create(out_dir_main, recursive = TRUE, showWarnings = FALSE)


out_dir_com <- "error_common/"
out_dir_gen <- "error_genewise/"

files_com <- list.files(out_dir_com, pattern = ".txt")
files_com

files_gen <- list.files(out_dir_gen, pattern = ".txt")
files_gen


pref_com <- strsplit2(gsub(pattern = "_est_common.txt", "", files_com), split = "_disp_")
pref_com
pref_gen <- strsplit2(gsub(pattern = "_est_genewise.txt", "", files_gen), split = "_disp_")
pref_gen

pref_sim <- intersect(pref_com[, 1], pref_gen[, 1])
pref_sim


##############################################################################
### Plots for common dispersion 
##############################################################################

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]



for(i in 1:length(pref_sim)){
  # i = 1
  
  dir_com <- paste0(out_dir_com, pref_sim[i], "_disp_", pref_com[1, 2], "_est_common.txt") 
  dir_gen <- paste0(out_dir_gen, pref_sim[i], "_disp_", pref_gen[1, 2], "_est_genewise.txt") 
  
  
  if(file.exists(dir_com) && file.exists(dir_gen)){
    
    out_dir <- paste0(out_dir_main, pref_sim[i], "_disp_combined")
    
    est_com <- read.table(dir_com, header = TRUE)
    est_gen <- read.table(dir_gen, header = TRUE)
    
    est_com$simulation <- "common"
    est_gen$simulation <- "genewise"
    
    est <- rbind(est_com, est_gen)
    est$simulation <- factor(est$simulation, levels = c("common", "genewise"))
    
    est$dispersion <- factor(est$dispersion, levels = c("genewise", "moderated", "common"))
    est$method <- factor(est$method, levels = c("ML-dirmult", "PL", "CR"))
    
    est$dispersionmethod <- interaction(est$dispersion, est$method)
    
    
    ### Absolute error
    
    error <- data.frame(error = est$est - est$true, simulation = est$simulation, dispersion = est$dispersion, method = est$method, dispersionmethod = est$dispersionmethod)
    
    error$error <- abs(error$error)
    
    ylim <- c(min(aggregate(. ~ dispersionmethod, error[, c("error", "dispersionmethod")], whisker_lower)[, 2]) - 1, max(aggregate(. ~ dispersionmethod, error[, c("error", "dispersionmethod")], whisker_upper)[, 2]) + 1)
    
    
    ggp <- ggplot(data = error, aes(y = error, x = dispersion, fill = method)) + 
      theme_bw() +
      ylab("Absolute error") +
      geom_boxplot(outlier.size = 0) +
      coord_cartesian(ylim = ylim) +
      theme(axis.text = element_text(size = 16), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
      facet_wrap(~ simulation)
    
    pdf(paste0(out_dir, "_boxplot_absolute.pdf"), 8, 5)
    print(ggp)
    dev.off()
    
    
    ### Error as a ratio
    
    error <- data.frame(error = est$est/est$true, simulation = est$simulation, dispersion = est$dispersion, method = est$method, dispersionmethod = est$dispersionmethod)
    
    
    ggp <- ggplot(data = error, aes(y = log10(error), x = dispersion, fill = method)) + 
      theme_bw() +
      ylab("Log10 of error ratio") +
      geom_boxplot(outlier.size = 1) +
      geom_hline(yintercept = 0, color="grey50", linetype = 2, size = 0.5) +
      theme(axis.text = element_text(size = 16), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
      facet_wrap(~ simulation)
    
    pdf(paste0(out_dir, "_boxplot_ratio_log.pdf"), 8, 5)
    print(ggp)
    dev.off()
    
    
  }
  
}


























