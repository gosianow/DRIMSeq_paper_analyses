######################################################
## ----- dispersion_fp_genewise_plots_run
## <<dispersion_fp_genewise_plots_run.R>>

# BioC 3.1
# Created 16 Nov 2015 

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

out_dir <- "fp_genewise/"

fp_files <- list.files(out_dir, pattern = ".txt")
fp_files



##############################################################################
### Plots 
##############################################################################


fp_files_list <- paste0(out_dir, fp_files)


for(i in 1:length(fp_files)){
  # i = 1
  
  out_dir <- file_path_sans_ext(fp_files_list[i])
  
  fp <- read.table(fp_files_list[i], header = TRUE)
  
  fp$dispersion <- factor(fp$dispersion, levels = c("genewise", "moderated", "common"))

  colnames(fp) <- c("fp", "Dispersion", "Method")

  ### fp
  
  ggp <- ggplot(data = fp, aes(y = fp, x = Dispersion)) + 
    theme_bw() +
    ylab("FP rate") +
    geom_boxplot(outlier.size = 1, fill = "grey80") +
    geom_hline(yintercept = 0.05, color="black", linetype = 2, size = 0.5) +
    theme(axis.text = element_text(size = 16), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16))
  
  pdf(paste0(out_dir, "_boxplot.pdf"), 5, 5)
  print(ggp)
  dev.off()
  

}


















































