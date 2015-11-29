######################################################
## ----- dispersion_error_moderation_real_plots_run
## <<dispersion_error_moderation_real_plots_run.R>>

# BioC 3.1
# Created 28 Nov 2015 

##############################################################################


library(ggplot2)
library(reshape2)
library(tools)
library(limma)
library(plyr)


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

out_dir_plots <- "error_moderation_real/"
dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)


##############################################################################
### Merge all results into one data frame
##############################################################################

sim_name='test2_'

n=c(3)


param_nm_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/nm_kim_kallisto_lognormal.txt'
param_nd_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/nd_common_kim_kallisto.txt'

param_pi_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/prop_kim_kallisto.txt'

param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/disp_genewise_kim_kallisto_lognormal.txt'


out_dir_res <- c("error_moderation_real/")

files <- list.files(out_dir_res[1], pattern = ".txt")
files


name_res <- c("est_moderation")




res_list <- list()
mse_list <- list()
ix <- 1

for(ix_n in 1:length(n)){
  # ix_n=1
  
  out_dir <- paste0(out_dir_res, sim_name, "n", n[ix_n], "_", basename(file_path_sans_ext(param_nm_path)), "_", basename(file_path_sans_ext(param_nd_path)), "_", basename(file_path_sans_ext(param_pi_path)), "_",  basename(file_path_sans_ext(param_gamma_path)), "_")
  
  out_name <- paste0(out_dir, name_res, ".txt")
  
  print(out_name)
  
  
  if(file.exists(out_name)){
    
    res_tmp <- read.table(out_name, header = TRUE, sep = "\t", as.is = TRUE)
    
    res_tmp$n <- n[ix_n]
    
    res_list[[ix]] <- res_tmp
    
    ### Splitting of results per run
    res_split <- split.data.frame(res_tmp, res_tmp$run)
    
    mse_tmp <- lapply(1:length(res_split), function(r){
      # r = 1
      rr <- res_split[[r]]
      rr$error_abs <- abs(rr$est - rr$true)
      rr$disp_prior_df <- factor(rr$disp_prior_df)
      
      out_mean <- aggregate(. ~ disp_prior_df, rr[, c("disp_prior_df", "error_abs")], mean)
      colnames(out_mean) <- c("disp_prior_df", "mean_error_abs")
      
      out_median <- aggregate(. ~ disp_prior_df, rr[, c("disp_prior_df", "error_abs")], median)
      colnames(out_median) <- c("disp_prior_df", "median_error_abs")
      
      out <- merge(out_mean, out_median, by = "disp_prior_df", sort = FALSE)
      
      out$n <- n[ix_n]
      
      return(out)
      
    })
    
    mse_list[[ix]] <- rbind.fill(mse_tmp)
    
    
    ix <- ix + 1
    
  }
  
}


res <- rbind.fill(res_list)

mse <- rbind.fill(mse_list)


##############################################################################
### Panel plots
##############################################################################

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]



### Adjust the order of the variables for plotting

res$disp_prior_df <- factor(res$disp_prior_df)



### Absolute error

error <- res[complete.cases(res), ]
error$error <- abs(error$est - error$true)


ggp <- ggplot(data = error, aes(y = error, x = disp_prior_df)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Absolute error") +
  xlab("Moderation") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) 

pdf(paste0(out_dir_plots, "error_absolute_violin.pdf"))
print(ggp)
dev.off()



ylim <- c(-2.5, 5)


ggp <- ggplot(data = error, aes(y = log10(error), x = disp_prior_df)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Log10 of absolute error") +
  xlab("Moderation") +
  # coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16))


pdf(paste0(out_dir_plots, "error_absolute_log_violin.pdf"))
print(ggp)
dev.off()




### Error as a ratio

error <- res[complete.cases(res), ]
error$error <- error$est/error$true

ylim <- c(-2, 2)

ggp <- ggplot(data = error, aes(y = log10(error), x = disp_prior_df)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5) +
  geom_hline(yintercept = 0, color="red", linetype = 2, size = 0.2) +
  theme_bw() +
  ylab("Log10 of error ratio") +
  xlab("Moderation") +
  # coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16))

pdf(paste0(out_dir_plots, "error_ratio_log_violin.pdf"))
print(ggp)
dev.off()







###### Plots of MSE (or mean absolute error)

### Adjust the order of the variables for plotting

mse$disp_prior_df <- factor(mse$disp_prior_df)


### plot mean 

ggp <- ggplot(data = mse, aes(y = mean_error_abs, x = disp_prior_df)) + 
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  # geom_hline(yintercept = log10(min_median), color="red", linetype = 2, size = 0.5) +
  theme_bw() +
  ylab("Mean absolute error") +
  xlab("Moderation") +
  # coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16))

pdf(paste0(out_dir_plots, "error_mean_absolute_boxplot.pdf"))
print(ggp)
dev.off()





### plot median 

ggp <- ggplot(data = mse, aes(y = median_error_abs, x = disp_prior_df)) + 
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  # geom_hline(yintercept = log10(min_median), color="red", linetype = 2, size = 0.5) +
  theme_bw() +
  ylab("Median absolute error") +
  xlab("Moderation") +
  # coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) 

pdf(paste0(out_dir_plots, "error_median_absolute_boxplot.pdf"))
print(ggp)
dev.off()












