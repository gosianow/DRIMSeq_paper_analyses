######################################################

# BioC 3.1
# Created 3 Nov 2015 

##############################################################################

library(pryr)
library(dirmult)
library(limma)
library(DRIMSeq)
library(ggplot2)

source("/home/gosia/R/drimseq_paper/simulations_dm/dm_simulate.R")

##############################################################################
# Read in the arguments
##############################################################################

rwd="/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/"
workers=2


# ## Read input arguments
# args <- (commandArgs(trailingOnly = TRUE))
# for (i in 1:length(args)) {
#   eval(parse(text = args[[i]]))
# }


print(rwd)
print(workers)



##############################################################################

dir.create(rwd, recursive = T, showWarnings = FALSE)
setwd(rwd)



##############################################################################
### Simulations for common dispersion MSE
##############################################################################

out_dir <- "./"


r <- 5 # Number of simulations
m <- 100 # Number of genes
n <- 3 # Number of samples
pi <- c(1/3, 1/3, 1/3) # Proportions
g0 <- 10 # Dispersion gamma0


mse_list <- lapply(1:r, function(i){
  # i = 1
  
  mse <- data.frame(matrix(0, nrow = 7, ncol = 3))
  colnames(mse) <- c("mse", "dispersion", "method")
  mse$dispersion <- c(rep(c("common", "genewise", "moderated"), 2), "genewise")
  mse$method <- c(rep(c("PL", "CR"), each = 3), "ML-dirmult")
  
  counts <- dm_simulate(m = m, n = n, pi = pi, g0 = g0, nm = 100, tot = "uni", nd = 0.25, mc.cores = workers)
  
  print(head(counts))
  
  group_split <- strsplit2(rownames(counts), ":")
  
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep("c1", ncol(counts)))
  
  disp_interval <- c(0, g0*100)
  
  ### No CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[1,"mse"] <- mean((common_dispersion(d) - g0)^2)
  
  mse[2,"mse"] <- mean((genewise_dispersion(d)$genewise_dispersion - g0)^2)
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = common_dispersion(d), disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[3,"mse"] <- mean((genewise_dispersion(d)$genewise_dispersion - g0)^2)
  
  
  ### With CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[4,"mse"] <- mean((common_dispersion(d) - g0)^2)
  
  mse[5,"mse"] <- mean((genewise_dispersion(d)$genewise_dispersion - g0)^2)
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = common_dispersion(d), disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[6,"mse"] <- mean((genewise_dispersion(d)$genewise_dispersion - g0)^2)
  
  
  ### gammas estimated with dirmult package
  
  dirmult_dispersion <- unlist(lapply(1:length(d@counts), function(g){
    # g = 1
    y <- t(d@counts[[g]])
    out <- dirmult(data = y, trace = FALSE, epsilon = 10^(4))
    sum(out$gamma)
    # (1-out$theta)/out$theta
    }))
  
  mse[7,"mse"] <- mean((dirmult_dispersion - g0)^2)
  
  return(mse)
  
})



mse <- do.call(rbind, mse_list)

write.table(mse, paste0(out_dir, "mse_common.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


################### Plots of MSE of dispersion gamma0



mse <- read.table(paste0(out_dir, "mse.txt"), header = TRUE)

mse$dispersion <- factor(mse$dispersion, levels = c("common", "genewise", "moderated"))
mse$method <- factor(mse$method, levels = c("ML-dirmult", "PL", "CR"))


ggp <- ggplot(data = mse, aes(y = log10(mse), x = dispersion, fill = method)) + 
  geom_boxplot()


pdf(paste0(out_dir, "mse_common_boxplot.pdf"))
print(ggp)
dev.off()




##############################################################################
### Simulations for common dispersion - error
##############################################################################

out_dir <- "./"


workers <- 5
r <- 100 # Number of simulations



m <- 100 # Number of genes
n <- 6 # Number of samples
pi <- c(1/3, 1/3, 1/3) # Proportions
g0 <- 10 # Dispersion gamma0


error_list <- lapply(1:r, function(i){
  # i = 1
  
  mse <- list()
  
  counts <- dm_simulate(m = m, n = n, pi = pi, g0 = g0, nm = 100, tot = "uni", nd = 0.25, mc.cores = workers)
  print(head(counts))
  
  group_split <- strsplit2(rownames(counts), ":")
  
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep("c1", ncol(counts)))
  
  disp_interval <- c(0, g0*100)
  
  ### No CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[1]] <- data.frame(mse = common_dispersion(d) - g0, dispersion = "common", method = "PL")
  
  mse[[2]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "genewise", method = "PL")
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = common_dispersion(d), disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[3]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "moderated", method = "PL")
  
  
  ### With CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[4]] <- data.frame(mse = common_dispersion(d) - g0, dispersion = "common", method = "CR")
  
  mse[[5]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "genewise", method = "CR")
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = common_dispersion(d), disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[6]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "moderated", method = "CR")
  
  
  ### gammas estimated with dirmult package
  
  dirmult_dispersion <- unlist(lapply(1:length(d@counts), function(g){
    # g = 1
    y <- t(d@counts[[g]])
    out <- dirmult(data = y, trace = FALSE, epsilon = 10^(4))
    sum(out$gamma)
    # (1-out$theta)/out$theta
  }))
  
  mse[[7]] <- data.frame(mse = dirmult_dispersion - g0, dispersion = "genewise", method = "ML-dirmult")
  
  mse <- do.call(rbind, mse)
  return(mse)
  
})




error <- do.call(rbind, arror_list)

write.table(error, paste0(out_dir, "error_common.txt"), quote = FALSE, sep = "\t", row.names = FALSE)


################### Plots of error of dispersion gamma0


error <- read.table(paste0(out_dir, "error.txt"), header = TRUE)

error$dispersion <- factor(error$dispersion, levels = c("common", "genewise", "moderated"))
error$method <- factor(error$method, levels = c("ML-dirmult", "PL", "CR"))

colnames(error) <- c("Error", "Dispersion", "Method")

ggp <- ggplot(data = error, aes(y = Error, x = Dispersion, fill = Method)) + 
  theme_bw() +
  geom_boxplot(outlier.size = 0) +
  coord_cartesian(ylim = quantile(error$Error, c(0.05, 0.95))) +
  geom_hline(yintercept=c(0), color="grey") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))


pdf(paste0(out_dir, "error_common_boxplot.pdf"))
print(ggp)
dev.off()




ggp <- ggplot(data = error, aes(y = Error, x = Dispersion, fill = Method)) + 
  theme_bw() +
  geom_violin(trim = FALSE) +
  coord_cartesian(ylim = quantile(error$Error, c(0.01, 0.98))) +
  geom_hline(yintercept=c(0), color="grey") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))


pdf(paste0(out_dir, "error_common_violin.pdf"))
print(ggp)
dev.off()





##############################################################################
### Simulations for genewise dispersion - error
##############################################################################

out_dir <- "./"


workers <- 5
r <- 50 # Number of simulations


m <- 100 # Number of genes
n <- 6 # Number of samples
pi <- c(1/3, 1/3, 1/3) # Proportions
g0_shape <- 10
g0_scale <- 2


### Choose a distribution for gamma0
pdf(paste0(out_dir, "gamma_distribution.pdf"))
x <- seq(1:50)
plot(x, dgamma(x, shape = g0_shape, scale = g0_scale, log = FALSE), type = "l", ylab = "", main = paste0("Gamma density with shape = ", g0_shape, " and scale = ", g0_scale))
dev.off()





error_list <- lapply(1:r, function(i){
  # i = 1
  
  mse <- list()
  
  g0 <- rgamma(m, shape = 10, scale = 2)
  
  counts <- dm_simulate(m = m, n = n, pi = pi, g0 = g0, nm = 100, tot = "uni", nd = 0.25, mc.cores = workers)
  print(head(counts))
  
  group_split <- strsplit2(rownames(counts), ":")
  
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep("c1", ncol(counts)))
  
  disp_interval <- c(0, max(g0)*100)
  
  ### No CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[1]] <- data.frame(mse = common_dispersion(d) - g0, dispersion = "common", method = "PL")
  
  mse[[2]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "genewise", method = "PL")
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = common_dispersion(d), disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[3]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "moderated", method = "PL")
  
  
  ### With CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[4]] <- data.frame(mse = common_dispersion(d) - g0, dispersion = "common", method = "CR")
  
  mse[[5]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "genewise", method = "CR")
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = common_dispersion(d), disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[6]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, dispersion = "moderated", method = "CR")
  
  
  ### gammas estimated with dirmult package
  
  dirmult_dispersion <- unlist(lapply(1:length(d@counts), function(g){
    # g = 1
    y <- t(d@counts[[g]])
    out <- dirmult(data = y, trace = FALSE, epsilon = 10^(4))
    sum(out$gamma)
    # (1-out$theta)/out$theta
  }))
  
  mse[[7]] <- data.frame(mse = dirmult_dispersion - g0, dispersion = "genewise", method = "ML-dirmult")
  
  mse <- do.call(rbind, mse)
  
  print(mem_used())
  rm("d")
  gc()
  print(mem_used())
  
  return(mse)
  
})




error <- do.call(rbind, error_list)

write.table(error, paste0(out_dir, "error_genewise.txt"), quote = FALSE, sep = "\t", row.names = FALSE)



################### Plots of error of dispersion gamma0

library(ggplot2)
whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]


error <- read.table(paste0(out_dir, "error_genewise.txt"), header = TRUE)

error$dispersion <- factor(error$dispersion, levels = c("common", "genewise", "moderated"))
error$method <- factor(error$method, levels = c("ML-dirmult", "PL", "CR"))

colnames(error) <- c("Error", "Dispersion", "Method")

ggp <- ggplot(data = error, aes(y = Error, x = Dispersion, fill = Method)) + 
  theme_bw() +
  geom_boxplot(outlier.size = 0) +
  coord_cartesian(ylim = c(min(aggregate(. ~ interaction(error$Dispersion, error$Method), error, whisker_lower)[, 2]) - 1, max(aggregate(. ~ interaction(error$Dispersion, error$Method), error, whisker_upper)[, 2]) + 1)) +
  geom_hline(yintercept=c(0), color="grey") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))


pdf(paste0(out_dir, "error_genewise_boxplot.pdf"))
print(ggp)
dev.off()




ggp <- ggplot(data = error, aes(y = Error, x = Dispersion, fill = Method)) + 
  theme_bw() +
  geom_violin(trim = FALSE) +
  # coord_cartesian(ylim = quantile(error$Error, c(0, 0.98))) +
  geom_hline(yintercept=c(0), color="grey") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))

pdf(paste0(out_dir, "error_genewise_violin.pdf"))
print(ggp)
dev.off()



##############################################################################
### Simulations for moderated dispersion 
##############################################################################

out_dir <- "./"


workers <- 5
r <- 5 # Number of simulations



m <- 10 # Number of genes
n <- 6 # Number of samples
pi <- c(1/3, 1/3, 1/3) # Proportions


### Choose a distribution for gamma0

# Gamma
g0_shape <- 10
g0_scale <- 2

pdf(paste0(out_dir, "gamma_distribution.pdf"))
x <- seq(1:50)
plot(x, dgamma(x, shape = g0_shape, scale = g0_scale, log = FALSE), type = "l", ylab = "", main = paste0("Gamma density with shape = ", g0_shape, " and scale = ", g0_scale))
dev.off()



# # Log-normal
# g0_location <- log(10) # median = exp(location)
# g0_scale <- log(1.5)
# 
# pdf(paste0(out_dir, "lognormal_distribution.pdf"))
# x <- seq(1:100)
# plot(x, dlnorm(x, meanlog = g0_location, sdlog = g0_scale, log = FALSE), type = "l", ylab = "", main = paste0("Log-normal density with location = ", round(g0_location, 2), " and scale = ", round(g0_scale, 2)))
# dev.off()


disp_prior_df <- seq(0, 1, by = 0.05)


error_list <- lapply(1:r, function(i){
  # i = 1
  
  g0 <- rgamma(m, shape = 10, scale = 2)
  
  counts <- dm_simulate(m = m, n = n, pi = pi, g0 = g0, nm = 100, tot = "uni", nd = 0.25, mc.cores = workers)
  print(head(counts))
  
  group_split <- strsplit2(rownames(counts), ":")
  
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep("c1", ncol(counts)))
  
  disp_interval <- c(0, max(g0)*100)

  ### With CR adjustement
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = FALSE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  common_dispersion <- common_dispersion(d)
  
  mse <- lapply(1:length(disp_prior_df), function(j){
    # j = 1
    
    d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = common_dispersion, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "common", disp_prior_df = disp_prior_df[j], disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
    
    out <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, disp_prior_df = disp_prior_df[j])
    
    })

  mse <- do.call(rbind, mse)
  return(mse)
  
})




error <- do.call(rbind, error_list)

write.table(error, paste0(out_dir, "error_moderation.txt"), quote = FALSE, sep = "\t", row.names = FALSE)



################### Plots of error of dispersion gamma0

library(ggplot2)

error <- read.table(paste0(out_dir, "error_moderation.txt"), header = TRUE)

error$disp_prior_df <- factor(error$disp_prior_df)
colnames(error) <- c("Error", "G0")

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]

ggp <- ggplot(data = error, aes(y = Error, x = G0)) + 
  theme_bw() +
  geom_boxplot(outlier.size = 0, fill = "grey80") +
  coord_cartesian(ylim = c(min(aggregate(. ~ G0, error, whisker_lower)[, 2]) - 1, max(aggregate(. ~ G0, error, whisker_upper)[, 2]) + 1)) +
  geom_hline(yintercept=c(0), color="grey") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))


pdf(paste0(out_dir, "error_moderation_boxplot.pdf"))
print(ggp)
dev.off()




ggp <- ggplot(data = error, aes(y = Error, x = G0)) + 
  theme_bw() +
  geom_violin(trim = FALSE, fill = "grey80") +
  # coord_cartesian(ylim = quantile(error$Error, c(0, 0.98))) +
  geom_hline(yintercept=c(0), color="grey") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "bottom", legend.title = element_text(size = 16), legend.text = element_text(size = 16))

pdf(paste0(out_dir, "error_moderation_violin.pdf"))
print(ggp)
dev.off()






##############################################################################
### Simulations for genewise dispersion - different optimization methods 
##############################################################################

out_dir <- "./"


workers <- 5
r <- 5 # Number of simulations



m <- 10 # Number of genes
n <- 6 # Number of samples
pi <- c(1/3, 1/3, 1/3) # Proportions


### Choose a distribution for gamma0

# Gamma
g0_shape <- 10
g0_scale <- 2

pdf(paste0(out_dir, "gamma_distribution.pdf"))
x <- seq(1:50)
plot(x, dgamma(x, shape = g0_shape, scale = g0_scale, log = FALSE), type = "l", ylab = "", main = paste0("Gamma density with shape = ", g0_shape, " and scale = ", g0_scale))
dev.off()




error_list <- lapply(1:r, function(i){
  # i = 1
  
  g0 <- rgamma(m, shape = 10, scale = 2)
  
  counts <- dm_simulate(m = m, n = n, pi = pi, g0 = g0, nm = 100, tot = "uni", nd = 0.25, mc.cores = workers)
  print(head(counts))
  
  group_split <- strsplit2(rownames(counts), ":")
  
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = rep("c1", ncol(counts)))
  
  disp_interval <- c(0, max(g0)*100)
  
  mse <- list()
  
  ### With CR adjustement
  
  ## grid
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = TRUE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "grid", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))

  mse[[1]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, method = "grid")

  ## optimize
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "optimize", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[2]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, method = "optimize")
  
  ## optim
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "optim", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[3]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, method = "optim")
  
  ## constrOptim
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = TRUE, disp_mode = "constrOptim", disp_interval = disp_interval, disp_tol = 1e-08, disp_init = 100, disp_init_weirMoM = TRUE, disp_grid_length = 21, disp_grid_range = c(-10, 10), disp_moderation = "none", disp_prior_df = 1, disp_span = 0.3, prop_mode = "constrOptimG", prop_tol = 1e-12, verbose = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = workers))
  
  mse[[4]] <- data.frame(mse = genewise_dispersion(d)$genewise_dispersion - g0, method = "constrOptim")
  
  mse <- do.call(rbind, mse)
  return(mse)
  
})




error <- do.call(rbind, error_list)

write.table(error, paste0(out_dir, "error_optimization.txt"), quote = FALSE, sep = "\t", row.names = FALSE)



################### Plots of error of dispersion gamma0

library(ggplot2)

error <- read.table(paste0(out_dir, "error_optimization.txt"), header = TRUE)

error$method <- factor(error$method, levels = c("grid", "optimize", "optim", "constrOptim"))
colnames(error) <- c("Error", "Method")

whisker_upper <- function(x) boxplot.stats(x)$stats[5]
whisker_lower <- function(x) boxplot.stats(x)$stats[1]

ggp <- ggplot(data = error, aes(y = Error, x = Method, fill = Method)) + 
  theme_bw() +
  geom_boxplot(outlier.size = 0) +
  coord_cartesian(ylim = c(min(aggregate(. ~ Method, error, whisker_lower)[, 2]) - 1, max(aggregate(. ~ Method, error, whisker_upper)[, 2]) + 1)) +
  geom_hline(yintercept=c(0), color="grey") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "none", legend.title = element_text(size = 16), legend.text = element_text(size = 16))


pdf(paste0(out_dir, "error_optimization_boxplot.pdf"))
print(ggp)
dev.off()




ggp <- ggplot(data = error, aes(y = Error, x = Method, fill = Method)) + 
  theme_bw() +
  geom_violin(trim = FALSE) +
  # coord_cartesian(ylim = quantile(error$Error, c(0, 0.98))) +
  geom_hline(yintercept=c(0), color="grey") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.position = "none", legend.title = element_text(size = 16), legend.text = element_text(size = 16))

pdf(paste0(out_dir, "error_optimization_violin.pdf"))
print(ggp)
dev.off()

















