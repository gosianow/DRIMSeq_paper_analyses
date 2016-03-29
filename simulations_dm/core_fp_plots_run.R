######################################################
## <<core_fp_plots_run.R>>

# BioC 3.2
# Created 8 Mar 2016

##############################################################################
Sys.time()
##############################################################################


library(ggplot2)
library(reshape2)
library(tools)
library(limma)
library(plyr)

##############################################################################
# Arguments for testing the code
##############################################################################

# rwd='/home/gosia/multinomial_project/simulations_dm/drimseq/'
# sim_name=''
# n=c(3,6,12)
# nm=c(100,1000)
# nd=0
# prop=c('prop_q3_uniform','prop_q3_kim_kallisto_fcutoff','prop_q10_uniform','prop_q10_kim_kallisto_fcutoff')
# disp=c('disp_common_kim_kallisto','disp_genewise_kim_kallisto_lognormal')
# pdf_width=15
# pdf_height=10
# fig_name='all_'


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


##############################################################################

setwd(rwd)

out_dir_res <- "core_fp/run/"
out_dir_plots <- "core_fp/"

out_suffix <- "core_fp"


##############################################################################
### Merge all results into one data frame
##############################################################################



res_list <- list()
mse_list <- list()
fp_list <- list()
ix <- 1

for(ix_n in 1:length(n)){
  
  for(ix_nm in 1:length(nm)){
    
    for(ix_prop in 1:length(prop)){
      
      for(ix_disp in 1:length(disp)){
        # ix_n=1; ix_nm=1; ix_prop=1; ix_disp=1
        
        out_name <- paste0(sim_name, "n", n[ix_n], "_nm", nm[ix_nm], "_nd", nd, "_", prop[ix_prop], "_", disp[ix_disp], "_")
        
        files <- list.files(out_dir_res, pattern = paste0(out_name, "est"))
        #print(files)
        
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
            rr$dispersion <- factor(rr$dispersion)
            rr$method <- factor(rr$method)
            
            out_mean <- aggregate(. ~ dispersion + method, rr[, c("dispersion", "method", "error_abs")], mean)
            colnames(out_mean) <- c("dispersion", "method", "mean_error_abs")
            
            out_median <- aggregate(. ~ dispersion + method, rr[, c("dispersion", "method", "error_abs")], median)
            colnames(out_median) <- c("dispersion", "method", "median_error_abs")
            
            out <- merge(out_mean, out_median, by = c("dispersion", "method"), sort = FALSE)
            
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
        #print(files)
        
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

res <- rbind.fill(res_list)
mse <- rbind.fill(mse_list)
fp <- rbind.fill(fp_list)


##############################################################################
### Panel plots
##############################################################################

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]



### Adjust the order of the variables for plotting

res$simulation <- "common"
res$simulation[grepl("genewise", res$disp)] <- "genewise"

res$simulation <- factor(res$simulation, levels = c("common", "genewise"))

res$dispersion <- factor(res$dispersion, levels = c("genewise", "moderated", "common"))
res$method <- factor(res$method, levels = c("true", "ML-dirmult", "PL", "CR"))

res$prop <- factor(res$prop, levels = prop)
levels(res$prop)

res$n <- factor(res$n, levels = n, labels = paste0("n=", n))
res$nm <- factor(res$nm, levels = nm, labels = paste0("m=", nm))

res$n_nm <- interaction(res$n, res$nm, lex.order = TRUE)
levels(res$n_nm)

res$n_nm_simulation <- interaction(res$n_nm, res$simulation, lex.order = TRUE)
levels(res$n_nm_simulation)

res$all_interactions <- interaction(res$dispersion, res$method, res$prop, drop = TRUE)





### Absolute error

error <- res[res$method != "true", ]
error$error <- abs(error$est - error$true)


ylim <- c(min(aggregate(. ~ all_interactions, error[, c("error", "all_interactions")], whisker_lower)[, "error"]) - 1, max(aggregate(. ~ all_interactions, error[, c("error", "all_interactions")], whisker_upper)[, "error"]) + 1)


ggp <- ggplot(data = error, aes(y = error, x = dispersion, fill = method)) + 
  geom_boxplot(outlier.size = 0.3, outlier.colour = "grey") +
  theme_bw() +
  ylab("Absolute error") +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(prop ~ n_nm_simulation)

pdf(paste0(out_dir_plots, "/", fig_name, "error_absolute_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()

png(paste0(out_dir_plots, "/", fig_name, "error_absolute_boxplot.png"), width = 100*pdf_width, height = 100*pdf_height)
print(ggp)
dev.off()



# ### Error as a ratio
# 
# error <- res
# error$error <- res$est/res$true
# 
# 
# ggp <- ggplot(data = error, aes(y = log10(error), x = dispersion, fill = method)) + 
#   geom_boxplot(outlier.size = 0.3, outlier.colour = "red") +
#   theme_bw() +
#   ylab("Log10 of error ratio") +
#   coord_cartesian(ylim = c(-2, 2.5)) +
#   geom_hline(yintercept = 0, color="grey70", linetype = 1, size = 0.3) +
#   theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16), panel.grid.major.x =element_blank()) +
#   geom_vline(xintercept = c(1.5, 2.5), color = "grey90", size = 0.05) +
#   facet_grid(prop ~ n_nm_simulation)
# 
# pdf(paste0(out_dir_plots, "/", fig_name, "error_ratio_log_boxplot.pdf"), width = pdf_width, height = pdf_height)
# print(ggp)
# dev.off()



### Estimates

## Plot true estimates with grey color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

ggn <- sum(!grepl("true", levels(res$method)))



ggp <- ggplot(data = res, aes(y = log10(est), x = dispersion, fill = method)) + 
  geom_violin(trim = TRUE, scale = "width", position = position_dodge(width = 0.9), width = 0.9) +
  geom_boxplot(outlier.size = NA, outlier.colour = NULL, position = position_dodge(width = 0.9), width = 0.5) +
  theme_bw() +
  ylab("Log10 of gamma_+") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  scale_fill_manual(values = c("grey", gg_color_hue(ggn))) +
  facet_grid(prop ~ n_nm_simulation)


pdf(paste0(out_dir_plots, "/", fig_name, "est_log_violin.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()



###### Plots of MSE (or mean absolute error)

### Adjust the order of the variables for plotting

mse <- mse[mse$method != "true", ]

mse$simulation <- "common"
mse$simulation[grepl("genewise", mse$disp)] <- "genewise"

mse$simulation <- factor(mse$simulation, levels = c("common", "genewise"))

mse$dispersion <- factor(mse$dispersion, levels = c("genewise", "moderated", "common"))
mse$method <- factor(mse$method, levels = c("ML-dirmult", "PL", "CR"))

mse$prop <- factor(mse$prop, levels = prop)
levels(mse$prop)

mse$n <- factor(mse$n, levels = n, labels = paste0("n=", n))
mse$nm <- factor(mse$nm, levels = nm, labels = paste0("m=", nm))

mse$n_nm <- interaction(mse$n, mse$nm, lex.order = TRUE)
levels(mse$n_nm)

mse$n_nm_simulation <- interaction(mse$n_nm, mse$simulation, lex.order = TRUE)
levels(mse$n_nm_simulation)

mse$all_interactions <- interaction(mse$dispersion, mse$method, mse$prop, mse$disp, mse$n_nm, drop = TRUE)


### plot mean 

ylim <- c(min(aggregate(. ~ all_interactions, mse[, c("mean_error_abs", "all_interactions")], whisker_lower)[, "mean_error_abs"]) - 1, max(aggregate(. ~ all_interactions, mse[, c("mean_error_abs", "all_interactions")], whisker_upper)[, "mean_error_abs"]) + 1)


ggp <- ggplot(data = mse, aes(y = mean_error_abs, x = dispersion, colour = method)) + 
  geom_boxplot(outlier.size = 1) +
  theme_bw() +
  ylab("Mean absolute error") +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(prop ~ n_nm_simulation)

pdf(paste0(out_dir_plots, "/", fig_name, "error_mean_absolute_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()



### plot median 

ylim <- c(min(aggregate(. ~ all_interactions, mse[, c("median_error_abs", "all_interactions")], whisker_lower)[, "median_error_abs"]) - 1, max(aggregate(. ~ all_interactions, mse[, c("median_error_abs", "all_interactions")], whisker_upper)[, "median_error_abs"]) + 1)


ggp <- ggplot(data = mse, aes(y = median_error_abs, x = dispersion, colour = method)) + 
  geom_boxplot(outlier.size = 1) +
  theme_bw() +
  ylab("Median absolute error") +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(prop ~ n_nm_simulation)

pdf(paste0(out_dir_plots, "/", fig_name, "error_median_absolute_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()




### False positives

### Adjust the order of the variables for plotting

fp$simulation <- "common"
fp$simulation[grepl("genewise", fp$disp)] <- "genewise"

fp$simulation <- factor(fp$simulation, levels = c("common", "genewise"))

fp$dispersion <- factor(fp$dispersion, levels = c("genewise", "moderated", "common"))

fp$method <- factor(fp$method, levels = c("true", "ML-dirmult", "PL", "CR"))

fp$prop <- factor(fp$prop, levels = prop)
levels(fp$prop)

fp$n <- factor(fp$n, levels = n, labels = paste0("n=", n))
fp$nm <- factor(fp$nm, levels = nm, labels = paste0("m=", nm))

fp$n_nm <- interaction(fp$n, fp$nm, lex.order = TRUE)
levels(fp$n_nm)

fp$n_nm_simulation <- interaction(fp$n_nm, fp$simulation, lex.order = TRUE)
levels(fp$n_nm_simulation)

fp$all_interactions <- interaction(fp$dispersion, fp$prop, fp$disp, fp$n_nm, drop = TRUE)



## Plot true estimates with grey color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

ggn <- sum(!grepl("true", levels(fp$method)))



### plot

ylim <- c(0, max(fp$fp, na.rm = TRUE) + 0.01)


ggp <- ggplot(data = fp, aes(y = fp, x = dispersion, fill = method)) + 
  geom_boxplot(outlier.size = 1) +
  geom_hline(yintercept = 0.05, color="black", linetype = 2, size = 0.3) +
  theme_bw() +
  ylab("FP rate") +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  scale_fill_manual(values = c("grey", gg_color_hue(ggn))) +
  facet_grid(prop ~ n_nm_simulation)

pdf(paste0(out_dir_plots, "/", fig_name, "fp_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()






### Plot distributions of p-values


## Plot true estimates with grey color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

ggn <- sum(!grepl("true", levels(fp$method)))



ggp <- ggplot(data = res, aes(x = pvalue, linetype = dispersion, colour = method)) + 
  geom_density(alpha = 0.75, trim = TRUE, adjust = 0.5) +
  theme_bw() +
  ylab("Density") +
  xlab("P-values") +
  coord_cartesian(xlim = c(0, 1)) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16)) +
  scale_colour_manual(values = c("grey", gg_color_hue(ggn))) +
  facet_grid(prop ~ n_nm_simulation)

pdf(paste0(out_dir_plots, "/", fig_name, "pvalues_density.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()





ggp <- ggplot(data = res, aes(x = pvalue, linetype = dispersion, colour = method)) + 
  geom_freqpoly(alpha = 0.75, breaks = seq(0, 1, 0.05)) +
  theme_bw() +
  ylab("Count") +
  xlab("P-values") +
  coord_cartesian(xlim = c(0, 1)) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16)) +
  scale_colour_manual(values = c("grey", gg_color_hue(ggn))) +
  facet_grid(prop ~ n_nm_simulation)

pdf(paste0(out_dir_plots, "/", fig_name, "pvalues_freqpoly.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()










































