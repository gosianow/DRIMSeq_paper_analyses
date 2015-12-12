######################################################
## ----- dispersion_error_filtering_plots_run
## <<dispersion_error_filtering_plots_run.R>>

# BioC 3.1
# Created 1 Dec 2015 

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
sim_name='run'
n=c(3)
nm=c(10000)
nd=0
param_pi_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/kim_kallisto/prop_q15_kim_kallisto_overall.txt'
param_gamma_path='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/kim_kallisto/disp_common_kim_kallisto.txt'
out_name_plots=''

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

out_dir_plots <- "error_filtering/"
dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)


##############################################################################
### Merge all results into one data frame
##############################################################################

out_dir_res <- c("error_filtering/run/")

suffix_res_est <- c("est_filtering")
suffix_res_fp <- c("fp_filtering")


prop <- basename(file_path_sans_ext(param_pi_path))
disp <- basename(file_path_sans_ext(param_gamma_path))



res_list <- list()
mse_list <- list()
fp_list <- list()
ix <- 1


for(ix_n in 1:length(n)){
  
  for(ix_nm in 1:length(nm)){
    
    for(ix_prop in 1:length(prop)){
      
      # ix_n=1; ix_nm=1; ix_prop=1;
      
      out_name <- paste0("n", n[ix_n], "_nm", nm[ix_nm], "_nd", nd, "_", basename(file_path_sans_ext(param_pi_path[ix_prop])), "_",  basename(file_path_sans_ext(param_gamma_path)), "_")
      
      ### est
      
      files <- list.files(out_dir_res, pattern = paste0(out_name, suffix_res_est, ".txt"))
      print(files)
      
      if(length(files) > 0){
        
        res_tmp_list <- list()
        mse_tmp_list <- list()
        
        for(i in 1:length(files)){
          # i = 1
          
          rr <- read.table(paste0(out_dir_res, files[i]), header = TRUE, sep = "\t", as.is = TRUE)
          rr$run <- i
          res_tmp_list[[i]] <- rr
          
          rr$error_abs <- abs(rr$est - rr$true)
          rr$max_features <- factor(rr$max_features)
          
          out_mean <- aggregate(. ~ max_features, rr[, c("max_features", "error_abs")], mean)
          colnames(out_mean) <- c("max_features", "mean_error_abs")
          
          out_median <- aggregate(. ~ max_features, rr[, c("max_features", "error_abs")], median)
          colnames(out_median) <- c("max_features", "median_error_abs")
          
          out <- merge(out_mean, out_median, by = "max_features", sort = FALSE)
          
          mse_tmp_list[[i]] <- out
          
        }
        
        res_tmp <- rbind.fill(res_tmp_list)
        res_tmp$n <- n[ix_n]
        res_tmp$nm <- nm[ix_nm]
        res_tmp$prop <- prop[ix_prop]
        res_list[[ix]] <- res_tmp
        
        mse_tmp <- rbind.fill(mse_tmp_list)
        mse_tmp$n <- n[ix_n]
        mse_tmp$nm <- nm[ix_nm]
        mse_tmp$prop <- prop[ix_prop]
        mse_list[[ix]] <- mse_tmp
      }
      ### fp
      
      files <- list.files(out_dir_res, pattern = paste0(out_name, suffix_res_fp, ".txt"))
      print(files)
      
      if(length(files) > 0){
        
        fp_tmp_list <- list()
        
        for(i in 1:length(files)){
          # i = 1
          
          rr <- read.table(paste0(out_dir_res, files[i]), header = TRUE, sep = "\t", as.is = TRUE)
          fp_tmp_list[[i]] <- rr
          
          
        }
        
        fp_tmp <- rbind.fill(fp_tmp_list)
        fp_tmp$n <- n[ix_n]
        fp_tmp$nm <- nm[ix_nm]
        fp_tmp$prop <- prop[ix_prop]
        fp_list[[ix]] <- fp_tmp
        
      }
      
      ix <- ix + 1
      
    }
  }
}


res <- rbind.fill(res_list)
mse <- rbind.fill(mse_list)
fp <- rbind.fill(fp_list)

##############################################################################
### Panel plots
##############################################################################

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]



max_features_levels <- sort(unique(res$max_features), decreasing = TRUE)
max_features_levels


### Adjust the order of the variables for plotting

res$max_features <- factor(res$max_features, levels = max_features_levels)

res$prop <- factor(res$prop, levels = prop)
res$n <- factor(res$n, levels = n, labels = paste0("n", n))
res$nm <- factor(res$nm, levels = nm, labels = paste0("nm", nm))

res$n_nm <- interaction(res$n, res$nm, lex.order = TRUE)
levels(res$n_nm)


### Absolute error

error <- res[complete.cases(res), ]
error$error <- abs(error$est - error$true)


ylim <- c(-2.5, 5)


ggp <- ggplot(data = error, aes(y = log10(error), x = max_features)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Log10 of absolute error") +
  xlab("Max features") +
  # coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(prop ~ n_nm)


pdf(paste0(out_dir_plots, out_name_plots, "error_absolute_log_violin.pdf"))
print(ggp)
dev.off()


### Estimates

true_disp <- res$true[1]

ggp <- ggplot(data = res, aes(y = log10(est), x = max_features)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  geom_hline(yintercept = log10(true_disp), color="black", linetype = 2, size = 0.5) +
  theme_bw() +
  ylab("Log 10 of gamma_+") +
  xlab("Max features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(prop ~ n_nm)


pdf(paste0(out_dir_plots, out_name_plots, "est_log_violin.pdf"))
print(ggp)
dev.off()



###### Plots of MSE (or mean absolute error)

### Adjust the order of the variables for plotting

mse$max_features <- factor(mse$max_features, levels = max_features_levels)
levels(mse$max_features)

mse$prop <- factor(mse$prop, levels = prop)
mse$n <- factor(mse$n, levels = n, labels = paste0("n", n))
mse$nm <- factor(mse$nm, levels = nm, labels = paste0("nm", nm))

mse$n_nm <- interaction(mse$n, mse$nm, lex.order = TRUE)
levels(mse$n_nm)


### plot mean 

ggp <- ggplot(data = mse, aes(y = mean_error_abs, x = max_features)) + 
  geom_boxplot(outlier.size = 1, fill = "grey80", width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Mean absolute error") +
  xlab("Max features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(prop ~ n_nm)

pdf(paste0(out_dir_plots, out_name_plots, "error_mean_absolute_boxplot.pdf"))
print(ggp)
dev.off()





### plot median 

ggp <- ggplot(data = mse, aes(y = median_error_abs, x = max_features)) + 
  geom_boxplot(outlier.size = 1, fill = "grey80", width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Median absolute error") +
  xlab("Max features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(prop ~ n_nm)

pdf(paste0(out_dir_plots, out_name_plots, "error_median_absolute_boxplot.pdf"))
print(ggp)
dev.off()







### False positives


fp$max_features <- factor(fp$max_features, levels = max_features_levels)
levels(fp$max_features)

fp$prop <- factor(fp$prop, levels = prop)
fp$n <- factor(fp$n, levels = n, labels = paste0("n", n))
fp$nm <- factor(fp$nm, levels = nm, labels = paste0("nm", nm))

fp$n_nm <- interaction(fp$n, fp$nm, lex.order = TRUE)
levels(fp$n_nm)



ylim <- c(0, max(fp$fp, na.rm = TRUE) + 0.01)


ggp <- ggplot(data = fp, aes(y = fp, x = max_features)) + 
  geom_boxplot(outlier.size = 1, fill = "grey60") +
  geom_hline(yintercept = 0.05, color="black", linetype = 2, size = 0.3) +
  theme_bw() +
  ylab("FP rate") +
  xlab("Max features") +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(prop ~ n_nm)

pdf(paste0(out_dir_plots, out_name_plots, "fp_boxplot.pdf"))
print(ggp)
dev.off()




























