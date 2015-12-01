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
library(DRIMSeq)

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/'
n=3 # Number of samples
param_nm_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/kim_kallisto/nm_kim_kallisto_lognormal.txt'
### Common dispersion of gene expression
param_nd_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/kim_kallisto/nd_common_kim_kallisto.txt'
param_pi_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/kim_kallisto/prop_kim_kallisto.txt'
### Genewise dispersion of feature proportions
param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/kim_kallisto/disp_genewise_kim_kallisto_lognormal.txt'

##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(n)
print(param_nm_path)
print(param_nd_path)
print(param_pi_path)
print(param_gamma_path)


##############################################################################

setwd(rwd)

out_dir_plots <- "error_moderation_real/"
dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)


##############################################################################
### Merge all results into one data frame
##############################################################################

out_dir_res <- c("error_moderation_real/run/")

suffix_res <- c("est_moderation")


res_list <- list()
mse_list <- list()
ix <- 1

for(ix_n in 1:length(n)){
  # ix_n=1
  
  out_name <- paste0("n", n[ix_n], "_", basename(file_path_sans_ext(param_nm_path)), "_", basename(file_path_sans_ext(param_nd_path)), "_", basename(file_path_sans_ext(param_pi_path)), "_",  basename(file_path_sans_ext(param_gamma_path)), "_")
  
  files <- list.files(out_dir_res, pattern = paste0(out_name, suffix_res, ".txt"))
  print(files)
  
  if(length(files) > 0){
    
    res_tmp_list <- list()
    mse_tmp_list <- list()
    
    for(i in 1:length(files)){
      
      rr <- read.table(paste0(out_dir_res, files[i]), header = TRUE, sep = "\t", as.is = TRUE)
      rr$run <- i
      res_tmp_list[[i]] <- rr
      
      rr$error_abs <- abs(rr$est - rr$true)
      rr$disp_prior_df <- factor(sprintf("%g", rr$disp_prior_df))
      
      out_mean <- aggregate(. ~ disp_prior_df, rr[, c("disp_prior_df", "error_abs")], mean)
      colnames(out_mean) <- c("disp_prior_df", "mean_error_abs")
      
      out_median <- aggregate(. ~ disp_prior_df, rr[, c("disp_prior_df", "error_abs")], median)
      colnames(out_median) <- c("disp_prior_df", "median_error_abs")
      
      out <- merge(out_mean, out_median, by = "disp_prior_df", sort = FALSE)
      
      mse_tmp_list[[i]] <- out
      
    }
    
    res_tmp <- rbind.fill(res_tmp_list)
    res_tmp$n <- n[ix_n]
    res_list[[ix]] <- res_tmp
    
    mse_tmp <- rbind.fill(mse_tmp_list)
    mse_tmp$n <- n[ix_n]
    mse_list[[ix]] <- mse_tmp
    
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


out_name <- paste0(basename(file_path_sans_ext(param_nm_path)), "_", basename(file_path_sans_ext(param_nd_path)), "_", basename(file_path_sans_ext(param_pi_path)), "_",  basename(file_path_sans_ext(param_gamma_path)), "_")

### Adjust the order of the variables for plotting

disp_prior_df_levels <- sprintf("%g", sort(unique(res$disp_prior_df)))

res$disp_prior_df <- factor(sprintf("%g", res$disp_prior_df), levels = disp_prior_df_levels)

levels(res$disp_prior_df)



### Absolute error

error <- res[complete.cases(res), ]
error$error <- abs(error$est - error$true)


# ggp <- ggplot(data = error, aes(y = error, x = disp_prior_df)) + 
#   geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
#   geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
#   theme_bw() +
#   ylab("Absolute error") +
#   xlab("Moderation") +
#   theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) 
# 
# pdf(paste0(out_dir_plots, out_name, "error_absolute_violin.pdf"))
# print(ggp)
# dev.off()



ylim <- c(-2.5, 5)


ggp <- ggplot(data = error, aes(y = log10(error), x = disp_prior_df)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Log10 of absolute error") +
  xlab("Moderation") +
  # coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16))


pdf(paste0(out_dir_plots, out_name, "error_absolute_log_violin.pdf"))
print(ggp)
dev.off()




# ### Error as a ratio
# 
# error <- res[complete.cases(res), ]
# error$error <- error$est/error$true
# 
# ylim <- c(-2, 2)
# 
# ggp <- ggplot(data = error, aes(y = log10(error), x = disp_prior_df)) + 
#   geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
#   geom_boxplot(outlier.size = 1, fill = NA, width = 0.5) +
#   geom_hline(yintercept = 0, color="red", linetype = 2, size = 0.2) +
#   theme_bw() +
#   ylab("Log10 of error ratio") +
#   xlab("Moderation") +
#   # coord_cartesian(ylim = ylim) +
#   theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16))
# 
# pdf(paste0(out_dir_plots, out_name, "error_ratio_log_violin.pdf"))
# print(ggp)
# dev.off()







###### Plots of MSE (or mean absolute error)

### Adjust the order of the variables for plotting


mse$disp_prior_df <- factor(mse$disp_prior_df, levels = disp_prior_df_levels)


### plot mean 

ggp <- ggplot(data = mse, aes(y = mean_error_abs, x = disp_prior_df)) + 
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  # geom_hline(yintercept = log10(min_median), color="red", linetype = 2, size = 0.5) +
  theme_bw() +
  ylab("Mean absolute error") +
  xlab("Moderation") +
  # coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16))

pdf(paste0(out_dir_plots, out_name, "error_mean_absolute_boxplot.pdf"))
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

pdf(paste0(out_dir_plots, out_name, "error_median_absolute_boxplot.pdf"))
print(ggp)
dev.off()




#### Plot dispersion vesrus mean for different disp_prior_df

# take only a subset of res

res_sub <- res[res$run == 1, ]

res_true <- res_sub[res_sub$disp_prior_df == 0, , drop = FALSE]
res_true$disp_prior_df = "true"
res_true$est <- res_true$true

resgg <- rbind.fill(res_sub, res_true)

mean_expression <- resgg$nm
genewise_dispersion <- resgg$est
nr_features <- resgg$q
disp_prior_df <- factor(resgg$disp_prior_df, levels = c("true", disp_prior_df_levels))

nlevels(disp_prior_df)



df <- data.frame(mean_expression = log10(mean_expression + 1), dispersion = log10(genewise_dispersion), nr_features = nr_features, disp_prior_df = disp_prior_df)

df_quant <- min(quantile(na.omit(df$nr_features), probs = 0.95), 30)
breaks <- seq(2, df_quant, ceiling(df_quant/10))


ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion", colour = "nr_features" )) +
  theme_bw() +
  xlab("Log10 of mean expression") +
  ylab("Log10 of gamma_+") +
  geom_point(size = 1.5, alpha = 0.7, na.rm = TRUE, shape = 1) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14), legend.position = "top") +
  guides(colour = guide_colorbar(barwidth = 20, barheight = 0.5)) +
  scale_colour_gradient(limits = c(2, max(breaks)), breaks = breaks, low = "royalblue2", high="red2", name = "Number of features", na.value = "red2") +
  facet_wrap(~ disp_prior_df, nrow = 3)


pdf(paste0(out_dir_plots, out_name, "dispersion_versus_mean.pdf"), 20, 12)
print(ggp)
dev.off()

































