######################################################
## ----- filtering_sim_plots_run
## <<filtering_sim_plots_run.R>>

# BioC 3.2
# Created 1 Dec 2015 
# Modified 18 Dec 2015

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

# rwd='/home/gosia/multinomial_project/simulations_dm/drimseq'
# sim_name=''
# n=c(3,6)
# nm=c(10000,1000)
# nd=0
# disp='disp_common_kim_kallisto'
# out_suffix='proportions_decay'
# pdf_width=7
# pdf_height=7


##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(args)

print(rwd)
print(sim_name)
print(n)
print(nm)
print(nd)
print(disp)
print(out_suffix)

##############################################################################

setwd(rwd)

out_dir_res <- "proportions/run/"
out_dir_plots <- "proportions/"




##############################################################################
### Merge all results into one data frame
##############################################################################


res_list <- list()
mse_list <- list()
fp_list <- list()

ix <- 1

for(ix_n in 1:length(n)){
  
  for(ix_nm in 1:length(nm)){
    
    for(ix_disp in 1:length(disp)){
      # ix_n=1; ix_nm=1;  ix_disp=1
      
      out_name <- paste0(sim_name, "n", n[ix_n], "_nm", nm[ix_nm], "_nd", nd, "_", disp[ix_disp], "_")
      
      files <- list.files(out_dir_res, pattern = paste0(out_name, "est_", out_suffix))
      files
      
      if(length(files) > 0){
        
        res_tmp_list <- list()
        mse_tmp_list <- list()
        
        for(i in 1:length(files)){
          # i = 1
          rr <- read.table(paste0(out_dir_res, files[i]), header = TRUE, sep = "\t", as.is = TRUE)
          head(rr)
          
          rr$run <- i
          res_tmp_list[[i]] <- rr
          
          # calculate mse
          rr$error_abs <- abs(rr$est - rr$true)
          rr$nr_features <- factor(rr$nr_features)
          
          out_mean <- aggregate(. ~ nr_features, rr[, c("nr_features", "error_abs")], mean)
          colnames(out_mean) <- c( "nr_features", "mean_error_abs")
          
          out_median <- aggregate(. ~ nr_features, rr[, c("nr_features", "error_abs")], median)
          colnames(out_median) <- c("nr_features", "median_error_abs")
          
          out <- merge(out_mean, out_median, by = c("nr_features"), sort = FALSE)
          
          mse_tmp_list[[i]] <- out
          
        }
        
        res_tmp <- rbind.fill(res_tmp_list)
        res_tmp$n <- n[ix_n]
        res_tmp$nm <- nm[ix_nm]
        res_tmp$disp <- disp[ix_disp]
        res_list[[ix]] <- res_tmp
        
        mse_tmp <- rbind.fill(mse_tmp_list)
        mse_tmp$n <- n[ix_n]
        mse_tmp$nm <- nm[ix_nm]
        mse_tmp$disp <- disp[ix_disp]
        mse_list[[ix]] <- mse_tmp
        
      }
      
      files <- list.files(out_dir_res, pattern = paste0(out_name, "fp_", out_suffix))
      files
      
      if(length(files) > 0){
        
        fp_tmp_list <- list()
        
        for(i in 1:length(files)){
          # i = 1
          rr <- read.table(paste0(out_dir_res, files[i]), header = TRUE, sep = "\t", as.is = TRUE)
          
          rr$disp_estimator <- "disp_est"
          fp_tmp_list[[paste0("disp_est", i)]] <- rr
          
        }
        
        files <- list.files(out_dir_res, pattern = paste0(out_name, "fptruedisp_", out_suffix))
        files
        
        if(length(files) > 0){
          
          for(i in 1:length(files)){
            # i = 1
            rr <- read.table(paste0(out_dir_res, files[i]), header = TRUE, sep = "\t", as.is = TRUE)
            
            rr$disp_estimator <- "disp_true"
            fp_tmp_list[[paste0("disp_true", i)]] <- rr
            
          }
        }
        
        fp_tmp <- rbind.fill(fp_tmp_list)
        fp_tmp$n <- n[ix_n]
        fp_tmp$nm <- nm[ix_nm]
        fp_tmp$disp <- disp[ix_disp]
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



nr_features_levels <- sort(unique(res$nr_features), decreasing = FALSE)
nr_features_levels


### Adjust the order of the variables for plotting

res$nr_features <- factor(res$nr_features, levels = nr_features_levels)

res$n <- factor(res$n, levels = n, labels = paste0("n", n))
res$nm <- factor(res$nm, levels = nm, labels = paste0("nm", nm))

res$n_nm <- interaction(res$n, res$nm, lex.order = TRUE)
levels(res$n_nm)


### Absolute error

error <- res[complete.cases(res), ]
error$error <- abs(error$est - error$true)


ggp <- ggplot(data = error, aes(y = log10(error), x = nr_features)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Log10 of absolute error") +
  xlab("Number of features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)


pdf(paste0(out_dir_plots, out_suffix, "_error_absolute_log_violin.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()


### Estimates

true_disp <- res$true[1]

ggp <- ggplot(data = res, aes(y = log10(est), x = nr_features)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  geom_hline(yintercept = log10(true_disp), color="black", linetype = 2, size = 0.5) +
  theme_bw() +
  ylab("Log 10 of gamma_+") +
  xlab("Number of features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)


pdf(paste0(out_dir_plots, out_suffix, "_est_log_violin.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()



true_disp <- res$true[1]
res$fp <- factor(res$pvalue < 0.05, levels = c("TRUE", "FALSE"))

ggp <- ggplot(data = res, aes(y = log10(est), x = nr_features, fill = fp)) + 
  geom_violin(trim = FALSE, scale = "area", position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.size = 1, outlier.colour = NULL, width = 0.6, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = log10(true_disp), color="black", linetype = 2, size = 0.5) +
  theme_bw() +
  ylab("Log 10 of gamma_+") +
  xlab("Number of features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 14), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)


pdf(paste0(out_dir_plots, out_suffix, "_est_log_violin_fp.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()


###### Plots of MSE (or mean absolute error)

### Adjust the order of the variables for plotting

mse$nr_features <- factor(mse$nr_features, levels = nr_features_levels)
levels(mse$nr_features)

mse$n <- factor(mse$n, levels = n, labels = paste0("n", n))
mse$nm <- factor(mse$nm, levels = nm, labels = paste0("nm", nm))

mse$n_nm <- interaction(mse$n, mse$nm, lex.order = TRUE)
levels(mse$n_nm)


### plot mean 

ggp <- ggplot(data = mse, aes(y = mean_error_abs, x = nr_features)) + 
  geom_boxplot(outlier.size = 1, fill = "grey80", width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Mean absolute error") +
  xlab("Number of features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)

pdf(paste0(out_dir_plots, out_suffix, "_error_mean_absolute_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()





### plot median 

ggp <- ggplot(data = mse, aes(y = median_error_abs, x = nr_features)) + 
  geom_boxplot(outlier.size = 1, fill = "grey80", width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Median absolute error") +
  xlab("Number of features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)

pdf(paste0(out_dir_plots, out_suffix, "_error_median_absolute_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()







### False positives


fp$nr_features <- factor(fp$nr_features, levels = nr_features_levels)
levels(fp$nr_features)

fp$n <- factor(fp$n, levels = n, labels = paste0("n", n))
fp$nm <- factor(fp$nm, levels = nm, labels = paste0("nm", nm))

fp$n_nm <- interaction(fp$n, fp$nm, lex.order = TRUE)
levels(fp$n_nm)



ylim <- c(0, max(fp$fp, na.rm = TRUE) + 0.01)


ggp <- ggplot(data = fp, aes(y = fp, x = nr_features, fill = disp_estimator)) + 
  geom_boxplot(outlier.size = 1) +
  geom_hline(yintercept = 0.05, color="black", linetype = 2, size = 0.3) +
  theme_bw() +
  ylab("FP rate") +
  xlab("Number of features") +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)

pdf(paste0(out_dir_plots, out_suffix, "_fp_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()




























