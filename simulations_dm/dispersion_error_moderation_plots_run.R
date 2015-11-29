######################################################
## ----- dispersion_error_moderation_plots_run
## <<dispersion_error_moderation_plots_run.R>>

# BioC 3.1
# Created 12 Nov 2015 
# Modified 28 Nov 2015
# Add MSE (or mean absolute error) plots like in edgeR paper 

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

out_dir_plots <- "error_moderation/"
dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)


##############################################################################
### Merge all results into one data frame
##############################################################################

sim_name=''
nd=0

n=c(3,6)
nm=c(100,1000)

prop=c('prop_q3_uniform','prop_q3_kim_kallisto_overall','prop_q10_uniform','prop_q10_kim_kallisto_overall')

param_pi_path=paste0('/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/',prop,'.txt')


disp=c('disp_genewise_kim_kallisto_lognormal')

param_gamma_path=paste0('/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/dm_parameters/',disp,'.txt')

out_dir_res <- c("error_moderation/")

files <- list.files(out_dir_res[1], pattern = ".txt")
files

name_res <- c("est_moderation")
simulation <- c("genewise")




res_list <- list()
mse_list <- list()
ix <- 1

for(ix_n in 1:length(n)){
  
  for(ix_nm in 1:length(nm)){
    
    for(ix_prop in 1:length(param_pi_path)){
      
      for(ix_disp in 1:length(param_gamma_path)){
        # ix_n=1; ix_nm=1; ix_prop=1; ix_disp=1
        
        out_dir <- paste0(out_dir_res[ix_disp], sim_name, "n", n[ix_n], "_nm", nm[ix_nm], "_nd", nd, "_", basename(file_path_sans_ext(param_pi_path[ix_prop])), "_",  basename(file_path_sans_ext(param_gamma_path[ix_disp])), "_")
        
        out_name <- paste0(out_dir, name_res[ix_disp], ".txt")
        
        print(out_name)
        
        if(file.exists(out_name)){
          
          res_tmp <- read.table(out_name, header = TRUE, sep = "\t", as.is = TRUE)
          
          res_tmp$simulation <- simulation[ix_disp]
          res_tmp$proportions <- prop[ix_prop]
          res_tmp$nm <- nm[ix_nm]
          res_tmp$n <- n[ix_n]
          
          res_list[[ix]] <- res_tmp
          
          ### Manual splitting of results. There were 10 sim of 500 genes and 13 different df 
          res_split <- split.data.frame(res_tmp, ceiling(seq_len(nrow(res_tmp))/(13*500)))
          
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
            
            out$simulation <- simulation[ix_disp]
            out$proportions <- prop[ix_prop]
            out$nm <- nm[ix_nm]
            out$n <- n[ix_n]
            
            return(out)
            
          })
          
          mse_list[[ix]] <- rbind.fill(mse_tmp)
          
          
          ix <- ix + 1
          
        }
        
        
      }
    }
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


res$proportions <- factor(res$proportions, levels = prop)
levels(res$proportions)

res$n <- factor(res$n)
res$n <- factor(res$n, labels = paste0("n", levels(res$n)))

res$nm <- factor(res$nm)
res$nm <- factor(res$nm, labels = paste0("nm", levels(res$nm)))

res$n_nm <- interaction(res$n, res$nm, lex.order = TRUE)

levels(res$n_nm)



res$all_interactions <- interaction(res$n_nm, res$proportions, res$disp_prior_df , drop = TRUE)


### Absolute error

error <- res[complete.cases(res), ]
error$error <- abs(error$est - error$true)


# min_median <- min(aggregate(. ~ all_interactions, error[, c("all_interactions", "error")], median)[, "error"])


ylim <- c(min(aggregate(. ~ all_interactions, error[, c("all_interactions", "error")], whisker_lower)[, "error"]) - 1, max(aggregate(. ~ all_interactions, error[, c("all_interactions", "error")], whisker_upper)[, "error"]) + 1)



ggp <- ggplot(data = error, aes(y = error, x = disp_prior_df)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  # geom_hline(yintercept = log10(min_median), color="red", linetype = 2, size = 0.5) +
  theme_bw() +
  ylab("Absolute error") +
  xlab("Moderation") +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(proportions ~ n_nm)

pdf(paste0(out_dir_plots, "error_absolute_violin.pdf"), 15, 10)
print(ggp)
dev.off()



ylim <- c(-2.5, 5)


ggp <- ggplot(data = error, aes(y = log10(error), x = disp_prior_df)) + 
  geom_violin(trim = FALSE, fill = "grey80", colour = "grey80") +
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  # geom_hline(yintercept = log10(min_median), color="red", linetype = 2, size = 0.5) +
  theme_bw() +
  ylab("Log10 of absolute error") +
  xlab("Moderation") +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(proportions ~ n_nm)

pdf(paste0(out_dir_plots, "error_absolute_log_violin.pdf"), 15, 10)
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
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(proportions ~ n_nm)

pdf(paste0(out_dir_plots, "error_ratio_log_violin.pdf"), 15, 10)
print(ggp)
dev.off()







###### Plots of MSE (or mean absolute error)

### Adjust the order of the variables for plotting

mse$disp_prior_df <- factor(mse$disp_prior_df)


mse$proportions <- factor(mse$proportions, levels = prop)
levels(mse$proportions)

mse$n <- factor(mse$n)
mse$n <- factor(mse$n, labels = paste0("n", levels(mse$n)))

mse$nm <- factor(mse$nm)
mse$nm <- factor(mse$nm, labels = paste0("nm", levels(mse$nm)))

mse$n_nm <- interaction(mse$n, mse$nm, lex.order = TRUE)

levels(mse$n_nm)


### plot mean 

ggp <- ggplot(data = mse, aes(y = mean_error_abs, x = disp_prior_df)) + 
  geom_boxplot(outlier.size = 1, fill = NA, width = 0.5, outlier.colour = NULL) +
  # geom_hline(yintercept = log10(min_median), color="red", linetype = 2, size = 0.5) +
  theme_bw() +
  ylab("Mean absolute error") +
  xlab("Moderation") +
  # coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(proportions ~ n_nm)

pdf(paste0(out_dir_plots, "error_mean_absolute_boxplot.pdf"), 15, 10)
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
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(proportions ~ n_nm)

pdf(paste0(out_dir_plots, "error_median_absolute_boxplot.pdf"), 15, 10)
print(ggp)
dev.off()












