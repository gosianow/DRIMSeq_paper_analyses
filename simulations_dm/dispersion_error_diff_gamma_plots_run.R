######################################################
## ----- dispersion_error_diff_gamma_plots_run
## <<dispersion_error_diff_gamma_plots_run.R>>

# BioC 3.1
# Created 12 Dec 2015 

##############################################################################

library(ggplot2)
library(reshape2)
library(tools)
library(plyr)

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_dm/drimseq/'


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


out_dir_plots <- "drimseq_0_3_1/error_diff_gamma/"
dir.create(out_dir_plots, recursive = TRUE, showWarnings = FALSE)


##############################################################################
### Merge all results into one data frame
##############################################################################

sim_name=''
nd=0

n=c(3,6)
nm=c(100,1000)

prop=c('prop_q3_uniform','prop_q3_kim_kallisto_overall','prop_q10_uniform','prop_q10_kim_kallisto_overall')

param_pi_path=paste0('/home/gosia/multinomial_project/simulations_dm/drimseq/dm_parameters/',prop,'.txt')


out_dir_res <- c("drimseq_0_3_1/error_diff_gamma/")
name_res <- c("est_diff_gamma")
simulation <- c("common")


files <- list.files(out_dir_res, pattern = ".txt")
files




res_list <- list()
ix <- 1

for(ix_n in 1:length(n)){
  
  for(ix_nm in 1:length(nm)){
    
    for(ix_prop in 1:length(param_pi_path)){
      # ix_n=1; ix_nm=1; ix_prop=1
      
      out_dir <- paste0(out_dir_res, sim_name, "n", n[ix_n], "_nm", nm[ix_nm], "_nd", nd, "_", basename(file_path_sans_ext(param_pi_path[ix_prop])), "_")
      
      out_name <- paste0(out_dir, name_res, ".txt")
      out_name
      
      if(file.exists(out_name)){
        
        res_tmp <- read.table(out_name, header = TRUE, sep = "\t", as.is = TRUE)
        
        res_tmp$proportions <- prop[ix_prop]
        res_tmp$nm <- nm[ix_nm]
        res_tmp$n <- n[ix_n]
        
        res_list[[ix]] <- res_tmp
        
        ix <- ix + 1
        
      }
      
      
      
    }
  }
}


res <- rbind.fill(res_list)








##############################################################################
### Plots for genewise dispersion - different optimization methods 
##############################################################################

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]



### Adjust the order of the variables for plotting

res$true <- factor(res$true)


res$proportions <- factor(res$proportions, levels = prop)
levels(res$proportions)

res$n <- factor(res$n)
res$n <- factor(res$n, labels = paste0("n", levels(res$n)))

res$nm <- factor(res$nm)
res$nm <- factor(res$nm, labels = paste0("nm", levels(res$nm)))


res$n_nm <- interaction(res$n, res$nm, lex.order = TRUE)

levels(res$n_nm)




### Plot estimates


true <- unique(res[, 2:ncol(res)])

est <- res


ggp <- ggplot(data = est, aes(y = log10(est), x = true)) + 
  geom_violin(trim = FALSE, show_guide = FALSE, fill = "grey") +
  geom_point(data = true, aes(y = log10(as.numeric(as.character(true))), x = true), colour = "blue") +
  theme_bw() +
  ylab("Log10 of gamma_+") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(proportions ~ n_nm)


pdf(paste0(out_dir_plots, "estimates_log_violin.pdf"), 10, 10)
print(ggp)
dev.off()


### Plot error


error <- res
error$error <- error$est - as.numeric(as.character(error$true))

ggp <- ggplot(data = error, aes(y = log10(abs(error)), x = true)) + 
  geom_violin(trim = FALSE, show_guide = FALSE, fill = "grey") +
  theme_bw() +
  ylab("Log10 of absolute error") +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_blank(), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(proportions ~ n_nm)


pdf(paste0(out_dir_plots, "error_log_violin.pdf"), 10, 10)
print(ggp)
dev.off()
































