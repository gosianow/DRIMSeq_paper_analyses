##############################################################################

# BioC 3.1
# Created 6 Nov 2015 

##############################################################################

library(DRIMSeq)
library(ggplot2)
library(MASS)

rwd="/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/"
setwd(rwd)

out_dir <- "dm_parameters/"
dir.create(out_dir, recursive = T, showWarnings = FALSE)

#######################################
# Test parameters
#######################################

# Proportions
pi <- c(1/3, 1/3, 1/3) 

# Genewise dispersion gamma_+ from gamma distribution
g0_shape <- 10
g0_scale <- 2


pdf(paste0(out_dir, "shape_", round(g0_shape), "_scale_", round(g0_scale), ".pdf"))
x <- seq(1, 50)
plot(x, dgamma(x, shape = g0_shape, scale = g0_scale, log = FALSE), type = "l", ylab = "", main = paste0("Gamma density with shape = ", g0_shape, " and scale = ", g0_scale))
dev.off()


save(pi, file = paste0(out_dir, "pi_", paste0(gsub("\\.", "", sprintf( "%.2f", pi)), collapse = "_"), ".RData"))
save(g0_shape, g0_scale, file = paste0(out_dir, "shape_", round(g0_shape), "_scale_", round(g0_scale), ".RData"))



#######################################
# Parameters based on Kim data
#######################################


main_data_dir <- "/home/Shared/data/seq/kim_adenocarcinoma/"

method_out <- "drimseq_0_3_1"
model <- "model_full"
count_method <- c("htseq", "kallisto")[1]

data_dir <- paste0(main_data_dir, method_out, "/",  model, "/", count_method, "/")

disp <- "genewise"
disp_mode <- "grid"
disp_moderation <- "none"


data_name <- paste0(data_dir, "/drimseq_", disp, "_", disp_mode, "_", disp_moderation, "_")

load(paste0(data_name, "d.Rdata"))

########################### common dispersion

common_disp <- common_dispersion(d)

write.table(round(common_disp), file = paste0(out_dir, "disp_common_kim_", count_method, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


########################### genewise dispersion

genewise_disp <- genewise_dispersion(d)
genewise_disp <- genewise_disp[complete.cases(genewise_disp),]


whisker_up <- ceiling(boxplot.stats(genewise_disp$genewise_dispersion)$stats[5])
quant_up <- quantile(genewise_disp$genewise_dispersion, 0.95, na.rm = TRUE)


pdf(paste0(out_dir, "disp_genewise_kim_", count_method, "_hist.pdf"))
ggp <- ggplot(data = genewise_disp, aes(x = genewise_dispersion)) + 
  geom_density() +
  # coord_cartesian(xlim = c(-1, quant_up)) 
  xlim(-1, quant_up)
print(ggp)
dev.off()



pdf(paste0(out_dir, "disp_genewise_kim_", count_method, "_hist_log.pdf"))
ggp <- ggplot(data = genewise_disp, aes(x = log(genewise_dispersion))) + 
  geom_density() 
print(ggp)
dev.off()


### Fit gamma

params <- fitdistr(x = genewise_disp$genewise_dispersion[genewise_disp$genewise_dispersion < quant_up], dgamma, list(shape = 0, scale = 100), lower = 0.001)
params[[1]]
g0_shape <- params[[1]][1]
g0_scale <- params[[1]][2]


pdf(paste0(out_dir, "disp_genewise_kim_", count_method, "_hist_gamma.pdf"))
ggp <- ggplot(data = genewise_disp, aes(x = genewise_dispersion)) + 
  geom_density() +
  xlim(-1, quant_up) +
  geom_line(data = data.frame(x = seq(1, round(quant_up)), y = dgamma(seq(1, round(quant_up)), shape = g0_shape, scale = g0_scale, log = FALSE)), aes(x = x, y = y), colour = "red")

print(ggp)
dev.off()



### Fit lognormal

params <- fitdistr(x = genewise_disp$genewise_dispersion[genewise_disp$genewise_dispersion < quant_up], "lognormal")
params[[1]]
g0_meanlog <- params[[1]][1]
g0_sdlog <- params[[1]][2]

pdf(paste0(out_dir, "disp_genewise_kim_", count_method, "_hist_lognormal.pdf"))
ggp <- ggplot(data = genewise_disp, aes(x = genewise_dispersion)) + 
  geom_density() +
  xlim(-1, whisker_up) +
  geom_line(data = data.frame(x = seq(1, round(whisker_up)), y = dlnorm(seq(1, round(whisker_up)), meanlog = g0_meanlog, sdlog = g0_sdlog)), aes(x = x, y = y), colour = "red")

print(ggp)
dev.off()


write.table(round(params[[1]], 2), file = paste0(out_dir, "disp_genewise_kim_", count_method, "_lognormal.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)


### Fit normal

params <- fitdistr(x = log(genewise_disp$genewise_dispersion[genewise_disp$genewise_dispersion < quant_up]), "normal")
params[[1]]
g0_mean <- params[[1]][1]
g0_sd <- params[[1]][2]


pdf(paste0(out_dir, "disp_genewise_kim_", count_method, "_hist_normal.pdf"))
ggp <- ggplot(data = genewise_disp, aes(x = log(genewise_dispersion))) + 
  geom_density() +
  geom_line(data = data.frame(x = seq(min(log(genewise_disp$genewise_dispersion)), max(log(genewise_disp$genewise_dispersion)), by = 0.01), y = dnorm(seq(min(log(genewise_disp$genewise_dispersion)), max(log(genewise_disp$genewise_dispersion)), by = 0.01), mean = g0_mean, sd = g0_sd)), aes(x = x, y = y), colour = "red")

print(ggp)
dev.off()



########################### 




write.table(c(1/3, 1/3, 1/3), file = paste0(out_dir, "prop_q3_uniform.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)









































