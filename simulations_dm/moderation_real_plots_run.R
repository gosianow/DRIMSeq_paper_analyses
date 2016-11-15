######################################################
## ----- moderation_real_plots_run
## <<moderation_real_plots_run.R>>

# BioC 3.2
# Created 28 Nov 2015 
# Updated 15 Nov 2016

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

rwd='/home/gosia/multinomial_project/simulations_dm/drimseq/'
out_dir_plots='moderation_real_v2'
out_dir_res='moderation_real_v2/run'
sim_name=''
n=c(3,6) # Number of samples
nm=c('nm_brooks_kallisto_lognormal','nm_brooks_htseq_lognormal','nm_kim_kallisto_lognormal','nm_kim_htseq_lognormal')
nd=c('nd_common_brooks_kallisto','nd_common_brooks_htseq','nd_common_kim_kallisto','nd_common_kim_htseq')
prop=c('prop_brooks_kallisto_fcutoff','prop_brooks_htseq_fcutoff','prop_kim_kallisto_fcutoff','prop_kim_htseq_fcutoff')
disp=c('disp_genewise_brooks_kallisto_lognormal','disp_genewise_brooks_htseq_lognormal','disp_genewise_kim_kallisto_lognormal','disp_genewise_kim_htseq_lognormal')
data_name=c('brooks','kim')
count_method=c('kallisto','htseq')
out_suffix='moderation_real'
pdf_width=14
pdf_height=7
fig_name='all_'
strip_text_size=16
text_size=16

##############################################################################
# Read in the arguments
##############################################################################

rm(list = ls())

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

cat(paste0(args, collapse = "\n"), fill = TRUE)


print(rwd)
print(sim_name)
print(n)
print(nm)
print(nd)
print(prop)
print(disp)
print(data_name)
print(count_method)
print(out_suffix)

##############################################################################

setwd(rwd)


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
              rr <- read.table(paste0(out_dir_res, "/", files[i]), header = TRUE, sep = "\t", as.is = TRUE)
              head(rr)
              
              rr$run <- i
              rr$mean_expression <- rr$nm
              res_tmp_list[[i]] <- rr
              
              # calculate mse
              rr$error <- rr$est - rr$true
              rr$error_abs <- abs(rr$est - rr$true)
              rr$disp_prior_df <- factor(rr$disp_prior_df)
              
              out_mean <- aggregate(. ~ disp_prior_df, rr[, c("disp_prior_df", "error_abs")], mean)
              colnames(out_mean) <- c( "disp_prior_df", "mean_error_abs")
              
              out_median <- aggregate(. ~ disp_prior_df, rr[, c("disp_prior_df", "error_abs")], median)
              colnames(out_median) <- c("disp_prior_df", "median_error_abs")
              
              out_median_raw <- aggregate(. ~ disp_prior_df, rr[, c("disp_prior_df", "error")], median)
              colnames(out_median_raw) <- c("disp_prior_df", "median_error")
              
              out <- merge(out_mean, out_median, by = c("disp_prior_df"), sort = FALSE)
              out <- merge(out, out_median_raw, by = c("disp_prior_df"), sort = FALSE)
              
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
              
              rr <- read.table(paste0(out_dir_res, "/",  files[i]), header = TRUE, sep = "\t", as.is = TRUE)
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

disp_prior_df_levels <- sprintf("%g", sort(unique(res$disp_prior_df)))
disp_prior_df_levels


res$disp_prior_df <- sprintf("%g", res$disp_prior_df)


### Adjust the order of the variables for plotting

res$disp_prior_df <- factor(res$disp_prior_df, levels = disp_prior_df_levels)

res$data_name <- NA 
for(i in 1:length(data_name))
  res$data_name[grepl(data_name[i], res$prop)] <- data_name[i]
res$data_name <- factor(res$data_name, levels = unique(data_name))

res$count_method <- NA 
for(i in 1:length(count_method))
  res$count_method[grepl(count_method[i], res$prop)] <- count_method[i]
res$count_method <- factor(res$count_method, levels = unique(count_method))


res$n <- factor(res$n, levels = n, labels = paste0("n=", n))

res$interaction <- interaction(res$data_name, res$n, lex.order = FALSE)
levels(res$interaction)




# --------------------------------------------------------------------
### Absolute error
# --------------------------------------------------------------------
error <- res[complete.cases(res), ]
error$error <- abs(error$est - error$true)


ggp <- ggplot(data = error, aes(y = log10(error), x = disp_prior_df)) + 
  geom_violin(trim = TRUE, scale = "width", fill = "cadetblue4", colour = "cadetblue4") +
  geom_boxplot(outlier.size = NA, alpha = 0, width = 0.5) +
  theme_bw() +
  ylab("Log10 of absolute error") +
  xlab("Moderation") +
  theme(axis.text = element_text(size = text_size), axis.text.x = element_text(size = text_size, angle = 90, vjust = 0.5, hjust = 1), axis.title.y = element_text(size = text_size, face = "bold"), axis.title.x = element_text(size = text_size, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = text_size), strip.text = element_text(size = strip_text_size)) +
  facet_grid(count_method ~ interaction)


pdf(paste0(out_dir_plots, "/", fig_name, "error_absolute_log_violin.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()


# --------------------------------------------------------------------
### Estimates
# --------------------------------------------------------------------

ggp <- ggplot(data = res, aes(y = log10(est), x = disp_prior_df)) + 
  geom_violin(trim = TRUE, scale = "width", fill = "cadetblue4", colour = "cadetblue4") +
  geom_boxplot(outlier.size = NA, alpha = 0, width = 0.5) +
  theme_bw() +
  ylab("Log10 of gamma_+") +
  xlab("Moderation") +
  theme(axis.text = element_text(size = text_size), axis.text.x = element_text(size = text_size, angle = 90, vjust = 0.5, hjust = 1), axis.title.y = element_text(size = text_size, face = "bold"), axis.title.x = element_text(size = text_size, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = text_size), strip.text = element_text(size = strip_text_size)) +
  facet_grid(count_method ~ interaction, scales = "free")


pdf(paste0(out_dir_plots, "/", fig_name, "est_log_violin.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()


# --------------------------------------------------------------------
# Plots of MSE (or mean absolute error)
# --------------------------------------------------------------------
### Adjust the order of the variables for plotting


mse$disp_prior_df <- sprintf("%g", as.numeric(as.character(mse$disp_prior_df)))

mse$disp_prior_df <- factor(mse$disp_prior_df, levels = disp_prior_df_levels)

mse$data_name <- NA 
for(i in 1:length(data_name))
  mse$data_name[grepl(data_name[i], mse$prop)] <- data_name[i]
mse$data_name <- factor(mse$data_name, levels = unique(data_name))

mse$count_method <- NA 
for(i in 1:length(count_method))
  mse$count_method[grepl(count_method[i], mse$prop)] <- count_method[i]
mse$count_method <- factor(mse$count_method, levels = unique(count_method))



mse$n <- factor(mse$n, levels = n, labels = paste0("n=", n))

mse$interaction <- interaction(mse$data_name, mse$n, lex.order = FALSE)
levels(mse$interaction)



### plot mean 

ggp <- ggplot(data = mse, aes(y = mean_error_abs, x = disp_prior_df)) + 
  geom_boxplot(outlier.size = 1, fill = "grey80", width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Mean absolute error") +
  xlab("Moderation") +
  theme(axis.text = element_text(size = text_size), axis.text.x = element_text(size = text_size, angle = 90, vjust = 0.5, hjust = 1), axis.title.y = element_text(size = text_size, face = "bold"), axis.title.x = element_text(size = text_size, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = text_size), strip.text = element_text(size = strip_text_size)) +
  facet_grid(count_method ~ interaction, scales = "free")

pdf(paste0(out_dir_plots, "/",  fig_name, "error_mean_absolute_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()


### plot median 

ggp <- ggplot(data = mse, aes(y = median_error_abs, x = disp_prior_df)) + 
  geom_boxplot(outlier.size = 1, fill = "grey80", width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Median absolute error") +
  xlab("Moderation") +
  theme(axis.text = element_text(size = text_size), axis.text.x = element_text(size = text_size, angle = 90, vjust = 0.5, hjust = 1), axis.title.y = element_text(size = text_size, face = "bold"), axis.title.x = element_text(size = text_size, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = text_size), strip.text = element_text(size = strip_text_size)) +
  facet_grid(count_method ~ interaction, scales = "free")

pdf(paste0(out_dir_plots, "/",  fig_name, "error_median_absolute_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()



ggp <- ggplot(data = mse, aes(y = log10(median_error_abs), x = disp_prior_df)) + 
  geom_boxplot(outlier.size = 1, fill = "grey80", width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Log10 of median absolute error") +
  xlab("Moderation") +
  theme(axis.text = element_text(size = text_size), axis.text.x = element_text(size = text_size, angle = 90, vjust = 0.5, hjust = 1), axis.title.y = element_text(size = text_size, face = "bold"), axis.title.x = element_text(size = text_size, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = text_size), strip.text = element_text(size = strip_text_size)) +
  facet_grid(count_method ~ interaction, scales = "free")

pdf(paste0(out_dir_plots, "/",  fig_name, "error_median_absolute_log10_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()



### plot median of raw error

ggp <- ggplot(data = mse, aes(y = median_error, x = disp_prior_df)) + 
  geom_hline(yintercept = 0, color="black", linetype = 2, size = 0.2) +
  geom_boxplot(outlier.size = 1, fill = "grey80", width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Median error") +
  xlab("Moderation") +
  theme(axis.text = element_text(size = text_size), axis.text.x = element_text(size = text_size, angle = 90, vjust = 0.5, hjust = 1), axis.title.y = element_text(size = text_size, face = "bold"), axis.title.x = element_text(size = text_size, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = text_size), strip.text = element_text(size = strip_text_size)) +
  facet_grid(count_method ~ interaction, scales = "free")

pdf(paste0(out_dir_plots, "/",  fig_name, "error_median_raw_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()





# --------------------------------------------------------------------
### False positives
# --------------------------------------------------------------------


fp$disp_prior_df <- sprintf("%g", as.numeric(as.character(fp$disp_prior_df)))

fp$disp_prior_df <- factor(fp$disp_prior_df, levels = disp_prior_df_levels)

fp$data_name <- NA 
for(i in 1:length(data_name))
  fp$data_name[grepl(data_name[i], fp$prop)] <- data_name[i]
fp$data_name <- factor(fp$data_name, levels = unique(data_name))

fp$count_method <- NA 
for(i in 1:length(count_method))
  fp$count_method[grepl(count_method[i], fp$prop)] <- count_method[i]
fp$count_method <- factor(fp$count_method, levels = unique(count_method))


fp$n <- factor(fp$n, levels = n, labels = paste0("n=", n))

fp$interaction <- interaction(fp$data_name, fp$n, lex.order = FALSE)
levels(fp$interaction)



ylim <- c(0, max(fp$fp, na.rm = TRUE) + 0.01)


ggp <- ggplot(data = fp, aes(y = fp, x = disp_prior_df)) + 
  geom_boxplot(outlier.size = 1, fill = "grey60") +
  geom_hline(yintercept = 0.05, color="black", linetype = 2, size = 0.3) +
  theme_bw() +
  ylab("FP rate") +
  xlab("Moderation") +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = text_size), axis.text.x = element_text(size = text_size, angle = 90, vjust = 0.5, hjust = 1), axis.title.y = element_text(size = text_size, face = "bold"), axis.title.x = element_text(size = text_size, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = text_size), strip.text = element_text(size = strip_text_size)) +
  facet_grid(count_method ~ interaction)

pdf(paste0(out_dir_plots, "/",  fig_name, "fp_boxplot.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()



# --------------------------------------------------------------------
### Plot distributions of p-values
# --------------------------------------------------------------------


ggp <- ggplot(data = res, aes(x = pvalue, colour = disp_prior_df)) +
  geom_density(alpha = 0.75, trim = TRUE, adjust = 0.75) +
  theme_bw() +
  ylab("Density") +
  xlab("P-values") +
  coord_cartesian(xlim = c(0, 1)) +
  theme(axis.text = element_text(size = text_size), axis.text.x = element_text(size = text_size), axis.title.y = element_text(size = text_size, face = "bold"), axis.title.x = element_text(size = text_size, face = "bold"), legend.position = "bottom", legend.title = element_text(size = text_size, face = "bold"), legend.text = element_text(size = text_size), strip.text = element_text(size = strip_text_size)) +
  scale_colour_manual(name = "Moderation", values = colorRampPalette(rev(c("red", "orange", "green", "darkgreen")))(nlevels(res$disp_prior_df)) ) +
  facet_grid(count_method ~ interaction) +
  guides(colour = guide_legend(nrow = 1)) 

pdf(paste0(out_dir_plots, "/", fig_name, "pvalues_density.pdf"), width = pdf_width, height = pdf_height)
print(ggp)
dev.off()





# --------------------------------------------------------------------
#### Plot dispersion vesrus mean for different disp_prior_df
# --------------------------------------------------------------------

# for(k in 1:length(n)){
#   
#   for(i in 1:length(data_name)){
#     
#     for(j in 1:length(count_method)){
#       # i = 1; j = 1
#       
#       
#       # take only a subset of res
#       
#       res_sub <- res[res$run == 1 & res$data_name == data_name[i] & res$count_method == count_method[j] & res$n == paste0("n", n[k]), ]
#       
#       res_true <- res_sub[res_sub$disp_prior_df == 0, , drop = FALSE]
#       res_true$disp_prior_df = "true"
#       res_true$est <- res_true$true
#       
#       resgg <- rbind.fill(res_sub, res_true)
#       
#       mean_expression <- resgg$mean_expression
#       genewise_dispersion <- resgg$est
#       nr_features <- resgg$q
#       disp_prior_df <- factor(resgg$disp_prior_df, levels = c("true", disp_prior_df_levels))
#       
#       nlevels(disp_prior_df)
#       
#       
#       
#       df <- data.frame(mean_expression = log10(mean_expression + 1), dispersion = log10(genewise_dispersion), nr_features = nr_features, disp_prior_df = disp_prior_df)
#       
#       df_quant <- min(quantile(na.omit(df$nr_features), probs = 0.95), 30)
#       breaks <- seq(2, df_quant, ceiling(df_quant/10))
#       
#       
#       ggp <- ggplot(df, aes_string(x = "mean_expression", y = "dispersion", colour = "nr_features" )) +
#         theme_bw() +
#         xlab("Log10 of mean expression") +
#         ylab("Log10 of gamma_+") +
#         geom_point(size = 1.5, alpha = 0.7, na.rm = TRUE, shape = 1) +
#         theme(axis.text = element_text(size=text_size), axis.title = element_text(size=18, face="bold"), legend.title = element_text(size=text_size, face="bold"), legend.text = element_text(size = text_size), legend.position = "top", strip.text = element_text(size = strip_text_size)) +
#         guides(colour = guide_colorbar(barwidth = 20, barheight = 0.5)) +
#         scale_colour_gradient(limits = c(2, max(breaks)), breaks = breaks, low = "royalblue2", high="red2", name = "Number of features", na.value = "red2") +
#         facet_wrap(~ disp_prior_df, nrow = 3)
#       
#       
#       pdf(paste0(out_dir_plots, "/","dispersion_versus_mean_n", n[k], "_", data_name[i], "_", count_method[j],".pdf"), width = pdf_width, height = pdf_height)
#       print(ggp)
#       dev.off()
#       
#       
#       
#       
#     }
#     
#   }
#   
# }




























#################################
# moderation_real_plots_run.R done!
#################################