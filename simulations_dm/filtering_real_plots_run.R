######################################################
## ----- filtering_real_plots_run
## <<filtering_real_plots_run.R>>

# BioC 3.2
# Created 16 Dec 2015 

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

rwd='/home/gosia/multinomial_project/simulations_dm/drimseq'
sim_name=''
n=3
nm=c('nm_brooks_kallisto_lognormal','nm_brooks_htseq_lognormal')
nd=c('nd_common_brooks_kallisto','nd_common_brooks_htseq')
prop=c('prop_brooks_kallisto_fcutoff','prop_brooks_htseq_fcutoff')
disp=c('disp_genewise_brooks_kallisto_lognormal','disp_genewise_brooks_htseq_lognormal')
data_name=c('brooks','brooks')
count_method=c('kallisto','htseq')

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
print(prop)
print(disp)
print(data_name)
print(count_method)


##############################################################################

setwd(rwd)

out_dir_res <- "filtering_real/run/"
out_dir_plots <- "filtering_real/"
out_suffix <- "filtering_real"


##############################################################################
### Merge all results into one data frame
##############################################################################


res_list <- list()
mse_list <- list()
fp_list <- list()
ix <- 1

for(ix_n in 1:length(n)){
  
  for(ix_nm in 1:length(nm)){
    
    for(ix_nd in 1:length(nd)){
      
      for(ix_prop in 1:length(prop)){
        
        for(ix_disp in 1:length(disp)){
          # ix_n=1; ix_nm=2; ix_nd=2; ix_prop=2; ix_disp=2
          
          out_name <- paste0(sim_name, "n",  n[ix_n], "_", nm[ix_nm], "_", nd[ix_nd], "_", prop[ix_prop], "_", disp[ix_disp], "_")
          
          files <- list.files(out_dir_res, pattern = paste0(out_name, "est"))
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
              rr$min_feature_expr <- factor(rr$min_feature_expr)
              
              out_mean <- aggregate(. ~ min_feature_expr, rr[, c("min_feature_expr", "error_abs")], mean)
              colnames(out_mean) <- c( "min_feature_expr", "mean_error_abs")
              
              out_median <- aggregate(. ~ min_feature_expr, rr[, c("min_feature_expr", "error_abs")], median)
              colnames(out_median) <- c("min_feature_expr", "median_error_abs")
              
              out <- merge(out_mean, out_median, by = c("min_feature_expr"), sort = FALSE)
              
              mse_tmp_list[[i]] <- out
              
            }
            
            res_tmp <- rbind.fill(res_tmp_list)
            res_tmp$n <- n[ix_n]
            res_tmp$nm <- nm[ix_nm]
            res_tmp$prop <- prop[ix_prop]
            res_tmp$disp <- disp[ix_disp]
            res_list[[ix]] <- res_tmp
            
            mse_tmp <- rbind.fill(mse_tmp_list)
            mse_tmp$n <- n[ix_n]
            mse_tmp$nm <- nm[ix_nm]
            mse_tmp$prop <- prop[ix_prop]
            mse_tmp$disp <- disp[ix_disp]
            mse_list[[ix]] <- mse_tmp
            
          }
          
          files <- list.files(out_dir_res, pattern = paste0(out_name, "fp"))
          files
          
          if(length(files) > 0){
            
            fp_tmp_list <- list()
            
            for(i in 1:length(files)){
              
              rr <- read.table(paste0(out_dir_res, files[i]), header = TRUE, sep = "\t", as.is = TRUE)
              fp_tmp_list[[i]] <- rr
              
              
            }
            
            fp_tmp <- rbind.fill(fp_tmp_list)
            fp_tmp$n <- n[ix_n]
            fp_tmp$nm <- nm[ix_nm]
            fp_tmp$prop <- prop[ix_prop]
            fp_tmp$disp <- disp[ix_disp]
            fp_list[[ix]] <- fp_tmp
            
          }
          
          
          ix <- ix + 1
        }
      }
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



min_feature_expr_levels <- sort(unique(res$min_feature_expr), decreasing = FALSE)
min_feature_expr_levels


### Adjust the order of the variables for plotting

res$min_feature_expr <- factor(res$min_feature_expr, levels = min_feature_expr_levels)

res$data_name <- NA 
for(i in 1:length(data_name))
res$data_name[grepl(data_name[i], res$prop)] <- data_name[i]
res$data_name <- factor(res$data_name, levels = unique(data_name))

res$count_method <- NA 
for(i in 1:length(count_method))
  res$count_method[grepl(count_method[i], res$prop)] <- count_method[i]
res$count_method <- factor(res$count_method, levels = unique(count_method))


### Absolute error

error <- res[complete.cases(res), ]
error$error <- abs(error$est - error$true)


ggp <- ggplot(data = error, aes(y = log10(error), x = min_feature_expr)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Log10 of absolute error") +
  xlab("Max features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(count_method ~ data_name)


pdf(paste0(out_dir_plots, "error_absolute_log_violin.pdf"), 7,7)
print(ggp)
dev.off()


### Estimates


ggp <- ggplot(data = res, aes(y = log10(est), x = min_feature_expr)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Log 10 of gamma_+") +
  xlab("Max features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(count_method ~ data_name)


pdf(paste0(out_dir_plots, "est_log_violin.pdf"),7,7)
print(ggp)
dev.off()



###### Plots of MSE (or mean absolute error)

### Adjust the order of the variables for plotting

mse$min_feature_expr <- factor(mse$min_feature_expr, levels = min_feature_expr_levels)

mse$data_name <- NA 
for(i in 1:length(data_name))
  mse$data_name[grepl(data_name[i], mse$prop)] <- data_name[i]
mse$data_name <- factor(mse$data_name, levels = unique(data_name))

mse$count_method <- NA 
for(i in 1:length(count_method))
  mse$count_method[grepl(count_method[i], mse$prop)] <- count_method[i]
mse$count_method <- factor(mse$count_method, levels = unique(count_method))


### plot mean 

ggp <- ggplot(data = mse, aes(y = mean_error_abs, x = min_feature_expr)) + 
  geom_boxplot(outlier.size = 1, fill = "grey80", width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Mean absolute error") +
  xlab("Max features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(count_method ~ data_name, scales = "free")

pdf(paste0(out_dir_plots, "error_mean_absolute_boxplot.pdf"),7,7)
print(ggp)
dev.off()


### plot median 

ggp <- ggplot(data = mse, aes(y = median_error_abs, x = min_feature_expr)) + 
  geom_boxplot(outlier.size = 1, fill = "grey80", width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Median absolute error") +
  xlab("Max features") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(count_method ~ data_name, scales = "free")

pdf(paste0(out_dir_plots, "error_median_absolute_boxplot.pdf"),7,7)
print(ggp)
dev.off()







### False positives

fp$min_feature_expr <- factor(fp$min_feature_expr, levels = min_feature_expr_levels)

fp$data_name <- NA 
for(i in 1:length(data_name))
  fp$data_name[grepl(data_name[i], fp$prop)] <- data_name[i]
fp$data_name <- factor(fp$data_name, levels = unique(data_name))

fp$count_method <- NA 
for(i in 1:length(count_method))
  fp$count_method[grepl(count_method[i], fp$prop)] <- count_method[i]
fp$count_method <- factor(fp$count_method, levels = unique(count_method))



ylim <- c(0, max(fp$fp, na.rm = TRUE) + 0.01)


ggp <- ggplot(data = fp, aes(y = fp, x = min_feature_expr)) + 
  geom_boxplot(outlier.size = 1, fill = "grey60") +
  geom_hline(yintercept = 0.05, color="black", linetype = 2, size = 0.3) +
  theme_bw() +
  ylab("FP rate") +
  xlab("Max features") +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(count_method ~ data_name)

pdf(paste0(out_dir_plots, "fp_boxplot.pdf"),7,7)
print(ggp)
dev.off()




























