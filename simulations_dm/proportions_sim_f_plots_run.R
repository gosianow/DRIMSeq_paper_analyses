######################################################
## ----- proportions_sim_f_plots_run
## <<proportions_sim_f_plots_run.R>>

# Created 17 Feb 2015 

##############################################################################

Sys.time()

##############################################################################

library(ggplot2)
library(reshape2)
library(tools)
library(limma)
library(plyr)
library(DRIMSeq)
library(RColorBrewer)


##############################################################################
# Arguments for testing the code
##############################################################################

# rwd='/home/gosia/multinomial_project/simulations_dm/drimseq'
# sim_name=''
# n=c(3,6)
# nm=c(1000,100000)
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

out_dir_res <- "proportions_f/run/"
out_dir_plots <- "proportions_f/"

strat_main <- "nr_features"
strat_extra <- c("disp_estimator", "test")

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
          rr[, strat_main] <- factor(rr[, strat_main])
          for(s in 1:length(strat_extra))
            rr[, strat_extra[s]] <- factor(rr[, strat_extra[s]])
          
          
          out_mean <- aggregate(as.formula(paste0(". ~ ", strat_main, " + ", paste0(strat_extra, collapse = " + "))), rr[, c(strat_main, strat_extra, "error_abs")], mean)
          colnames(out_mean) <- c(strat_main, strat_extra, "mean_error_abs")
          
          out_median <- aggregate(as.formula(paste0(". ~ ", strat_main, " + ", paste0(strat_extra, collapse = " + "))), rr[, c(strat_main, strat_extra, "error_abs")], median)
          colnames(out_median) <- c(strat_main, strat_extra, "median_error_abs")
          
          out <- merge(out_mean, out_median, by = c(strat_main, strat_extra), sort = FALSE)
          
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



strat_main_levels <- sort(unique(res[, strat_main]), decreasing = TRUE)
strat_main_levels


### Adjust the order of the variables for plotting

res[, strat_main] <- factor(res[, strat_main], levels = strat_main_levels)
res$strat_main <- res[, strat_main]

for(s in 1:length(strat_extra))
  res[, strat_extra[s]] <- factor(res[, strat_extra[s]])

res$disp_estimator <- relevel(res$disp_estimator, ref = "true")


res$n <- factor(res$n, levels = n, labels = paste0("n=", n))
res$nm <- factor(res$nm, levels = nm, labels = paste0("nm=", nm))

res$n_nm <- interaction(res$n, res$nm, lex.order = TRUE)
levels(res$n_nm)


### Absolute error

error <- res[res$disp_estimator != "true" & res$test == "lr", ]
error$error <- abs(error$est - error$true)


ggp <- ggplot(data = error, aes(y = log10(error), x = strat_main, fill = disp_estimator)) + 
  geom_violin(trim = FALSE, colour = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.size = 0.4, width = 0.5, outlier.colour = NULL, position = position_dodge(width = 0.8)) +
  theme_bw() +
  ylab("Log10 of absolute error") +
  xlab(strat_main) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)


pdf(paste0(out_dir_plots, out_suffix, "_error_absolute_log_violin.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()



### Estimates

true_disp <- res$true[1]

ggp <- ggplot(data = error, aes(y = log10(est), x = strat_main, fill = disp_estimator)) + 
  geom_violin(trim = FALSE, colour = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.size = 0.4, width = 0.5, outlier.colour = NULL, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = log10(true_disp), color="black", linetype = 2, size = 0.5) +
  theme_bw() +
  ylab("Log 10 of gamma_+") +
  xlab(strat_main) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)


pdf(paste0(out_dir_plots, out_suffix, "_est_log_violin.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()




###### Plots of MSE (or mean absolute error)

### Adjust the order of the variables for plotting
mse <- mse[mse$disp_estimator != "true" & mse$test == "lr", ]

mse[, strat_main] <- factor(mse[, strat_main], levels = strat_main_levels)
mse$strat_main <- mse[, strat_main]

for(s in 1:length(strat_extra))
  mse[, strat_extra[s]] <- factor(mse[, strat_extra[s]])

mse$n <- factor(mse$n, levels = n, labels = paste0("n=", n))
mse$nm <- factor(mse$nm, levels = nm, labels = paste0("nm=", nm))

mse$n_nm <- interaction(mse$n, mse$nm, lex.order = TRUE)
levels(mse$n_nm)


### plot mean 

ggp <- ggplot(data = mse, aes(y = mean_error_abs, x = strat_main, fill = disp_estimator)) + 
  geom_boxplot(outlier.size = 1, width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Mean absolute error") +
  xlab(strat_main) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)

pdf(paste0(out_dir_plots, out_suffix, "_error_mean_absolute_boxplot.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()


### plot median 

ggp <- ggplot(data = mse, aes(y = median_error_abs, x = strat_main, fill = disp_estimator)) + 
  geom_boxplot(outlier.size = 1, width = 0.5, outlier.colour = NULL) +
  theme_bw() +
  ylab("Median absolute error") +
  xlab(strat_main) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)

pdf(paste0(out_dir_plots, out_suffix, "_error_median_absolute_boxplot.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()






### False positives

fp[, strat_main] <- factor(fp[, strat_main], levels = strat_main_levels)
fp$strat_main <- fp[, strat_main]

for(s in 1:length(strat_extra))
  fp[, strat_extra[s]] <- factor(fp[, strat_extra[s]])

fp$disp_estimator <- relevel(fp$disp_estimator, ref = "true")

fp$n <- factor(fp$n, levels = n, labels = paste0("n=", n))
fp$nm <- factor(fp$nm, levels = nm, labels = paste0("nm=", nm))

fp$n_nm <- interaction(fp$n, fp$nm, lex.order = TRUE)
levels(fp$n_nm)


fp$test_disp_estimator <- interaction(fp$test, fp$disp_estimator, lex.order = TRUE)
levels(fp$test_disp_estimator)


ylim <- c(0, max(fp$fp, na.rm = TRUE) + 0.01)


ggp <- ggplot(data = fp, aes(y = fp, x = strat_main, fill = test_disp_estimator)) + 
  geom_boxplot(outlier.size = 1) +
  geom_hline(yintercept = 0.05, color="black", linetype = 2, size = 0.3) +
  theme_bw() +
  ylab("FP rate") +
  xlab(strat_main) +
  coord_cartesian(ylim = ylim) +
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 16)) +
  facet_grid(nm ~ n)

pdf(paste0(out_dir_plots, out_suffix, "_fp_boxplot.pdf"),width = pdf_width, height = pdf_height)
print(ggp)
dev.off()




### p-values

for(i in levels(res$test)){
  # i = "lr"
  
  ggp <- ggplot(data = res[res$test == i, ], aes(x = pvalue, linetype = disp_estimator, colour = strat_main)) + 
    geom_density(alpha = 0.75, trim = TRUE, adjust = 0.5) +
    theme_bw() +
    ylab("Density") +
    xlab("P-values") +
    coord_cartesian(xlim = c(0, 1)) +
    theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16)) +
    scale_colour_manual(name = strat_main, values = colorRampPalette(c("red", 'orange', 'blue'))(nlevels(res[, strat_main])) ) +
    scale_linetype_discrete(name = "Dispersion") +
    facet_grid(nm ~ n)
  
  pdf(paste0(out_dir_plots, out_suffix, "_pvalues_", i, "_density.pdf"), width = pdf_width, height = pdf_height)
  print(ggp)
  dev.off()
  
  
  
  ggp <- ggplot(data = res[res$test == i, ], aes(x = pvalue, linetype = disp_estimator, colour = strat_main)) + 
    geom_freqpoly(binwidth = 0.05, alpha = 0.75) +
    theme_bw() +
    ylab("Count") +
    xlab("P-values") +
    coord_cartesian(xlim = c(0, 1)) +
    theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16)) +
    scale_colour_manual(name = strat_main, values = colorRampPalette(c("red", 'orange', 'blue'))(nlevels(res[, strat_main])) ) +
    scale_linetype_discrete(name = "Dispersion") +
    facet_grid(nm ~ n)
  
  pdf(paste0(out_dir_plots, out_suffix, "_pvalues_", i, "_freqpoly.pdf"), width = pdf_width, height = pdf_height)
  print(ggp)
  dev.off()
  
  
  
}



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


### QQ plots


n <- 3
nm <- 1000

res <- res[res$n == paste0("n=", n) & res$nm == paste0("nm=", nm), ]



for(i in 1:nlevels(res$test)){
  
  ggpl <- list()
  
  for(j in 1:nlevels(res$nr_features)){
    # i = 2; j = 8
    
    test <- levels(res$test)[i]
    nr_features <- levels(res$nr_features)[j]
    
    
    switch(test, 
           
           f = {
             
             res_tmp <- res[res$test == test & res$nr_features == nr_features, ]
             
             print(all.equal(res_tmp[, "df1"], as.numeric(as.character(res_tmp[ , "nr_features"])) - 1 ))
             
             ggpl[[j]] <- ggplot(data = res_tmp, aes(sample = test_statistic, colour = disp_estimator)) + 
               geom_point(stat = "qq", distribution = qf, dparams = list(df1 = as.numeric(nr_features) - 1, df2 = (n - 1) * (as.numeric(nr_features) - 1))) +
               geom_abline(intercept = 0, slope = 1, color="black", linetype = 2, size = 0.3) +
               theme_bw() +
               theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16)) +
               ggtitle(paste0("Q-Q plot for genes with ", nr_features, " features"))
             
             
           },
           
           lr = {
             
             res_tmp <- res[res$test == test & res$nr_features == nr_features, ]
             
             print(all.equal(res_tmp[, "df"], as.numeric(as.character(res_tmp[ , "nr_features"])) - 1 ))
             
             ggpl[[j]] <- ggplot(data = res_tmp, aes(sample = test_statistic, colour = disp_estimator)) + 
               geom_point(stat = "qq", distribution = qchisq, dparams = list(df = as.numeric(nr_features) - 1)) +
               geom_abline(intercept = 0, slope = 1, color="black", linetype = 2, size = 0.3) +
               theme_bw() +
               theme(axis.text = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 16)) +
               ggtitle(paste0("Q-Q plot for genes with ", nr_features, " features"))
             
           }
           
    )
    
  }
  
  
  png(paste0(out_dir_plots, out_suffix, "_qqplot_", test, ".png"), width = 300 * pdf_width, height = 150 * pdf_height)
  print(multiplot(plotlist = ggpl, cols = 4))
  dev.off()
  
}











