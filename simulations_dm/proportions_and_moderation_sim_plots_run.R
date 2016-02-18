######################################################
## ----- proportions_and_moderation_sim_plots_run
## <<proportions_and_moderation_sim_plots_run.R>>

# Created 14 Jan 2015 

##############################################################################

Sys.time()

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
# nm=c(1000,100000)
# nd=0
# disp='disp_genewise_kim_kallisto_lognormal'
# out_suffix='proportions_uniform'
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

out_dir_res <- "proportions_and_moderation/run/"
out_dir_plots <- "proportions_and_moderation/"

feature_filter <- "nr_features"

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
      
      #       files <- list.files(out_dir_res, pattern = paste0(out_name, "est_", out_suffix))
      #       files
      
      pattern <- gsub("\\+", "\\\\+", paste0(out_name, "est_", out_suffix))
      pattern
      
      files <- list.files(path = out_dir_res, pattern = pattern)
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
          rr[, feature_filter] <- factor(rr[, feature_filter])
          rr$disp_estimator <- factor(rr$disp_estimator)
          
          
          out_mean <- aggregate(as.formula(paste0(". ~ ", feature_filter, "+ disp_estimator")), rr[, c(feature_filter, "disp_estimator", "error_abs")], mean)
          colnames(out_mean) <- c(feature_filter, "disp_estimator", "mean_error_abs")
          
          out_median <- aggregate(as.formula(paste0(". ~ ", feature_filter, "+ disp_estimator")), rr[, c(feature_filter, "disp_estimator", "error_abs")], median)
          colnames(out_median) <- c(feature_filter, "disp_estimator", "median_error_abs")
          
          out <- merge(out_mean, out_median, by = c(feature_filter, "disp_estimator"), sort = FALSE)
          
          mse_tmp_list[[i]] <- out
          
        }
        
        res_tmp <- rbind.fill(res_tmp_list)
        res_tmp$n <- n[ix_n]
        res_tmp$nm <- nm[ix_nm]
        res_tmp$disp <- disp[ix_disp]
        res_list[[paste0(ix)]] <- res_tmp
        
        mse_tmp <- rbind.fill(mse_tmp_list)
        mse_tmp$n <- n[ix_n]
        mse_tmp$nm <- nm[ix_nm]
        mse_tmp$disp <- disp[ix_disp]
        mse_list[[paste0(ix)]] <- mse_tmp
        
      }
      
#       files <- list.files(out_dir_res, pattern = paste0(out_name, "fp_", out_suffix))
#       files
      
      pattern <- gsub("\\+", "\\\\+", paste0(out_name, "fp_", out_suffix))
      pattern
      
      files <- list.files(path = out_dir_res, pattern = pattern)
      files
      
      if(length(files) > 0){
        
        fp_tmp_list <- list()
        
        for(i in 1:length(files)){
          # i = 1
          rr <- read.table(paste0(out_dir_res, files[i]), header = TRUE, sep = "\t", as.is = TRUE)
          head(rr)
          
          fp_tmp_list[[i]] <- rr
          
        }
        
        
        fp_tmp <- rbind.fill(fp_tmp_list)
        fp_tmp$n <- n[ix_n]
        fp_tmp$nm <- nm[ix_nm]
        fp_tmp$disp <- disp[ix_disp]
        fp_list[[paste0(ix)]] <- fp_tmp
        
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



feature_filter_levels <- sort(unique(res[, feature_filter]), decreasing = TRUE)
feature_filter_levels


### Adjust the order of the variables for plotting

res[, feature_filter] <- factor(res[, feature_filter], levels = feature_filter_levels)
res$feature_filter <- res[, feature_filter]

res$disp_estimator <- factor(res$disp_estimator)
res$disp_estimator <- relevel(res$disp_estimator, ref = "true")


res$n <- factor(res$n, levels = n, labels = paste0("n=", n))
res$nm <- factor(res$nm, levels = nm, labels = paste0("nm=", nm))

res$n_nm <- interaction(res$n, res$nm, lex.order = TRUE)
levels(res$n_nm)


### Absolute error

error <- res[res$disp_estimator != "true", ]
error$error <- abs(error$est - error$true)


ggp <- ggplot(data = error, aes(y = log10(error), x = feature_filter, fill = disp_estimator)) + 
  geom_violin(trim = FALSE, colour = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.size = 0.4, width = 0.5, outlier.colour = NULL, position = position_dodge(width = 0.8)) +
  theme_bw() +
  ylab("Log10 of absolute error") +
  xlab(feature_filter) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)


pdf(paste0(out_dir_plots, out_suffix, "_error_absolute_log_violin.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()



### Estimates

ggp <- ggplot(data = res, aes(y = log10(est), x = feature_filter, fill = disp_estimator)) + 
  geom_violin(trim = FALSE, colour = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.size = 0.4, width = 0.5, outlier.colour = NULL, position = position_dodge(width = 0.8)) +
  theme_bw() +
  ylab("Log 10 of gamma_+") +
  xlab(feature_filter) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)


pdf(paste0(out_dir_plots, out_suffix, "_est_log_violin.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()




###### Plots of MSE (or mean absolute error)

### Adjust the order of the variables for plotting
mse <- mse[mse$disp_estimator != "true", ]

mse[, feature_filter] <- factor(mse[, feature_filter], levels = feature_filter_levels)
mse$feature_filter <- mse[, feature_filter]

mse$disp_estimator <- factor(mse$disp_estimator)

mse$n <- factor(mse$n, levels = n, labels = paste0("n=", n))
mse$nm <- factor(mse$nm, levels = nm, labels = paste0("nm=", nm))

mse$n_nm <- interaction(mse$n, mse$nm, lex.order = TRUE)
levels(mse$n_nm)


### plot mean 

ggp <- ggplot(data = mse, aes(y = mean_error_abs, x = feature_filter, fill = disp_estimator)) + 
  geom_boxplot(outlier.size = 1, width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Mean absolute error") +
  xlab(feature_filter) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)

pdf(paste0(out_dir_plots, out_suffix, "_error_mean_absolute_boxplot.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()


### plot median 

ggp <- ggplot(data = mse, aes(y = median_error_abs, x = feature_filter, fill = disp_estimator)) + 
  geom_boxplot(outlier.size = 1, width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Median absolute error") +
  xlab(feature_filter) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)

pdf(paste0(out_dir_plots, out_suffix, "_error_median_absolute_boxplot.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()






### False positives

fp[, feature_filter] <- factor(fp[, feature_filter], levels = feature_filter_levels)
fp$feature_filter <- fp[, feature_filter]

fp$disp_estimator <- factor(fp$disp_estimator)
fp$disp_estimator <- relevel(fp$disp_estimator, ref = "true")

fp$n <- factor(fp$n, levels = n, labels = paste0("n=", n))
fp$nm <- factor(fp$nm, levels = nm, labels = paste0("nm=", nm))

fp$n_nm <- interaction(fp$n, fp$nm, lex.order = TRUE)
levels(fp$n_nm)



ylim <- c(0, max(fp$fp, na.rm = TRUE) + 0.01)


ggp <- ggplot(data = fp, aes(y = fp, x = feature_filter, fill = disp_estimator)) + 
  geom_boxplot(outlier.size = 1) +
  geom_hline(yintercept = 0.05, color="black", linetype = 2, size = 0.3) +
  theme_bw() +
  ylab("FP rate") +
  xlab(feature_filter) +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)

pdf(paste0(out_dir_plots, out_suffix, "_fp_boxplot.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()




### p-values



ggp <- ggplot(data = res, aes(x = pvalue, linetype = disp_estimator, colour = feature_filter)) + 
  geom_density(alpha = 0.75, trim = TRUE, adjust = 0.5) +
  theme_bw() +
  ylab("Density") +
  xlab("P-values") +
  coord_cartesian(xlim = c(0, 1)) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16)) +
  scale_colour_discrete(name = feature_filter) +
  scale_linetype_discrete(name = "Dispersion") +
  facet_grid(nm ~ n)

pdf(paste0(out_dir_plots, out_suffix, "_pvalues_density.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()



ggp <- ggplot(data = res, aes(x = pvalue, linetype = disp_estimator, colour = feature_filter)) + 
  geom_freqpoly(binwidth = 0.05, alpha = 0.75) +
  theme_bw() +
  ylab("Count") +
  xlab("P-values") +
  coord_cartesian(xlim = c(0, 1)) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16)) +
  scale_colour_discrete(name = feature_filter) +
  scale_linetype_discrete(name = "Dispersion") +
  facet_grid(nm ~ n)

pdf(paste0(out_dir_plots, out_suffix, "_pvalues_freqpoly.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()

























