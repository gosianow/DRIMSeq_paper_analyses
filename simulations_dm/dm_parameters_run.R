######################################################
## ----- dm_parameters_run
## <<dm_parameters_run.R>>

# BioC 3.1
# Created 6 Nov 2015 

##############################################################################

library(DRIMSeq)
library(ggplot2)
library(MASS)
library(edgeR)
library(reshape2)
library(matrixStats)
library(plyr)

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1'

count_method=c('htseq','kallisto')[2]

# main_data_dir='/home/Shared/data/seq/kim_adenocarcinoma'
# data_name='kim'

main_data_dir='/home/Shared/data/seq/brooks_pasilla'
data_name='brooks'


##############################################################################
# Read in the arguments
##############################################################################


## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(count_method)
print(main_data_dir)
print(data_name)

##############################################################################

setwd(rwd)

out_dir <- paste0("dm_parameters/", data_name, "_", count_method, "/")
dir.create(out_dir, recursive = T, showWarnings = FALSE)


##############################################################################
# Load data for the parameter estimation
##############################################################################

method_out <- "drimseq_0_3_1"

model='model_full'

load_data_dir <- paste0(main_data_dir, "/", method_out, "/",  model, "/", count_method, "/")

disp <- "genewise"
disp_mode <- "grid"
disp_moderation <- "none"

load_data_name <- paste0(load_data_dir, "/drimseq_", disp, "_", disp_mode, "_", disp_moderation, "_")

stopifnot(file.exists(paste0(load_data_name, "d.Rdata")))

load(paste0(load_data_name, "d.Rdata"))

use_sample <- as.character(samples(d)$sample_id[1])

use_sample

#######################################
### common dispersion
#######################################

common_disp <- common_dispersion(d)

write.table(round(common_disp), file = paste0(out_dir, "disp_common_",data_name,"_", count_method, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


#######################################
### genewise dispersion
#######################################

genewise_disp <- genewise_dispersion(d)
genewise_disp <- genewise_disp[complete.cases(genewise_disp),]


whisker_up <- ceiling(boxplot.stats(genewise_disp$genewise_dispersion)$stats[5])
quant_up <- quantile(genewise_disp$genewise_dispersion, 0.95, na.rm = TRUE)


pdf(paste0(out_dir, "disp_genewise_",data_name,"_", count_method, "_hist.pdf"))
ggp <- ggplot(data = genewise_disp, aes(x = genewise_dispersion)) + 
  geom_density() +
  # coord_cartesian(xlim = c(-1, quant_up)) 
  xlim(-1, quant_up)
print(ggp)
dev.off()



pdf(paste0(out_dir, "disp_genewise_",data_name,"_", count_method, "_hist_log.pdf"))
ggp <- ggplot(data = genewise_disp, aes(x = log(genewise_dispersion))) + 
  geom_density() 
print(ggp)
dev.off()


### Fit gamma

params <- fitdistr(x = genewise_disp$genewise_dispersion[genewise_disp$genewise_dispersion < quant_up], dgamma, list(shape = 0, scale = 100), lower = 0.001)
params[[1]]
g0_shape <- params[[1]][1]
g0_scale <- params[[1]][2]


pdf(paste0(out_dir, "disp_genewise_",data_name,"_", count_method, "_hist_gamma.pdf"))
ggp <- ggplot(data = genewise_disp, aes(x = genewise_dispersion)) + 
  geom_density() +
  xlim(-1, quant_up) +
  geom_line(data = data.frame(x = seq(1, round(quant_up)), y = dgamma(seq(1, round(quant_up)), shape = g0_shape, scale = g0_scale, log = FALSE)), aes(x = x, y = y), colour = "red") +
  ggtitle(paste0("Gamma fit with shape = ", round(g0_shape, 2), " and scale = ", round(g0_scale, 2)))

print(ggp)
dev.off()



### Fit lognormal

params <- fitdistr(x = genewise_disp$genewise_dispersion[genewise_disp$genewise_dispersion < quant_up], "lognormal")
params[[1]]
g0_meanlog <- params[[1]][1]
g0_sdlog <- params[[1]][2]

pdf(paste0(out_dir, "disp_genewise_",data_name,"_", count_method, "_hist_lognormal.pdf"))
ggp <- ggplot(data = genewise_disp, aes(x = genewise_dispersion)) + 
  geom_density() +
  xlim(-1, whisker_up) +
  geom_line(data = data.frame(x = seq(1, round(whisker_up)), y = dlnorm(seq(1, round(whisker_up)), meanlog = g0_meanlog, sdlog = g0_sdlog)), aes(x = x, y = y), colour = "red") +
  ggtitle(paste0("Lognormal fit with meanlog = ", round(g0_meanlog, 2), " and sdlog = ", round(g0_sdlog, 2)))

print(ggp)
dev.off()


write.table(round(params[[1]], 2), file = paste0(out_dir, "disp_genewise_",data_name,"_", count_method, "_lognormal.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)


### Fit normal

params <- fitdistr(x = log(genewise_disp$genewise_dispersion[genewise_disp$genewise_dispersion < quant_up]), "normal")
params[[1]]
g0_mean <- params[[1]][1]
g0_sd <- params[[1]][2]


pdf(paste0(out_dir, "disp_genewise_",data_name,"_", count_method, "_hist_normal.pdf"))
ggp <- ggplot(data = genewise_disp, aes(x = log(genewise_dispersion))) + 
  geom_density() +
  geom_line(data = data.frame(x = seq(min(log(genewise_disp$genewise_dispersion)), max(log(genewise_disp$genewise_dispersion)), by = 0.01), y = dnorm(seq(min(log(genewise_disp$genewise_dispersion)), max(log(genewise_disp$genewise_dispersion)), by = 0.01), mean = g0_mean, sd = g0_sd)), aes(x = x, y = y), colour = "red") +
  ggtitle(paste0("Normal fit with mean = ", round(g0_mean, 2), " and sd = ", round(g0_sd, 2)))

print(ggp)
dev.off()



#######################################
### dispersion for gene expression NB()
#######################################


metadata <- read.table(paste0(main_data_dir, "/3_metadata/metadata.xls"), stringsAsFactors = FALSE, sep="\t", header=TRUE) 
metadata_org <- metadata


count_dir <- paste0(main_data_dir, "/2_counts/", count_method, "/")

### load counts
counts_list <- lapply(1:length(metadata_org$sampleName), function(i){
  # i = 1
  cts <- read.table(paste0(count_dir, metadata_org$sampleName[i], ".counts"), header = FALSE, as.is = TRUE)
  colnames(cts) <- c("group_id", metadata_org$sampleName[i])  
  return(cts)
})

counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)
counts <- counts[!grepl(pattern = "_", counts$group_id),]


### Prepare data
group_split <- strsplit2(counts[,1], ":")
counts <- counts[, -1]
### order the samples like in metadata!!!
counts <- counts[, metadata_org$sampleName]


d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sampleName, group = metadata$condition)


gene_counts <- lapply(1:length(d@counts), function(g){ colSums(d@counts[[g]]) })
gene_counts <- do.call(rbind, gene_counts)
rownames(gene_counts) <- names(d)


dge <- DGEList(gene_counts, group = samples(d)$group)

dge <- calcNormFactors(dge)


dge$samples$lib.size


cpm <- cpm(dge, normalized.lib.sizes=TRUE)

keep <- rowSums(cpm > 1) >= ncol(dge)
dge <- dge[ keep, ]


pdf(paste0(out_dir, "nm_",data_name,"_", count_method, "_MDS.pdf"))
plotMDS(dge, col = c("darkgreen", "blue")[samples(d)$group])
dev.off()


dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

pdf(paste0(out_dir, "nm_",data_name,"_", count_method, "_BCV.pdf"))
plotBCV(dge)
dev.off()

dge$common.dispersion


write.table(round(dge$common.dispersion, 2), file = paste0(out_dir, "nd_common_",data_name,"_", count_method, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


bcv <- sqrt(dge$common.dispersion)
bcv



pdf(paste0(out_dir, "nd_tagwise_",data_name,"_", count_method, "_hist.pdf"))
ggp <- ggplot(data = data.frame(tagwise_dispersion = dge$tagwise.dispersion), aes(x = tagwise_dispersion)) + 
  geom_density() 
print(ggp)
dev.off()




##############################################################################
### gene expression
##############################################################################

gene_expr <- dge$counts[, use_sample]


whisker_up <- ceiling(boxplot.stats(gene_expr)$stats[5])
quant_up <- quantile(gene_expr, 0.95, na.rm = TRUE)



pdf(paste0(out_dir, "nm_",data_name,"_", count_method, "_hist.pdf"))
ggp <- ggplot(data.frame(gene_expr = gene_expr), aes(x = gene_expr)) +
  geom_density() +
  xlim(0, whisker_up)
print(ggp)
dev.off()


pdf(paste0(out_dir, "nm_",data_name,"_", count_method, "_hist_log.pdf"))
ggp <- ggplot(data.frame(gene_expr = gene_expr), aes(x = log(gene_expr))) +
  geom_density()
print(ggp)
dev.off()


### Fit negative binomial distribution to gene expression

params <- fitdistr(x = gene_expr[gene_expr <= whisker_up], "negative binomial")
params[[1]]

nm_size <- params[[1]][1]
nm_mu <- params[[1]][2]

ggp <- ggplot(data = data.frame(gene_expr = gene_expr), aes(x = gene_expr)) + 
  geom_density() +
  geom_line(data = data.frame(x = seq(1, round(whisker_up)), y = dnbinom(seq(1, round(whisker_up)), size = nm_size, mu = nm_mu)), aes(x = x, y = y), colour = "red") +
  ggtitle(paste0("Negative binomial fit with mu = ", round(nm_mu, 2), " and size = ", round(nm_size, 2))) +
  xlim(0, whisker_up)

pdf(paste0(out_dir, "nm_",data_name,"_", count_method, "_hist_negative_binomial.pdf"))
print(ggp)
dev.off()


write.table(round(params[[1]], 2), file = paste0(out_dir, "nm_",data_name,"_", count_method, "_negative_binomial.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)


### Fit lognormal distribution to gene expression

params <- fitdistr(x = gene_expr, "lognormal")
params[[1]]

nm_meanlog <- params[[1]][1]
nm_sdlog <- params[[1]][2]

ggp <- ggplot(data = data.frame(gene_expr = gene_expr), aes(x = gene_expr)) + 
  geom_density() +
  geom_line(data = data.frame(x = seq(1, round(whisker_up)), y = dlnorm(seq(1, round(whisker_up)), meanlog = nm_meanlog, sdlog = nm_sdlog)), aes(x = x, y = y), colour = "red") +
  ggtitle(paste0("Lognormal fit with meanlog = ", round(nm_meanlog, 2), " and sdlog = ", round(nm_sdlog, 2))) +
  xlim(0, whisker_up)

pdf(paste0(out_dir, "nm_",data_name,"_", count_method, "_hist_lognormal.pdf"))
print(ggp)
dev.off()


write.table(round(params[[1]], 2), file = paste0(out_dir, "nm_",data_name,"_", count_method, "_lognormal.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)


### Fit normal distribution to log of gene expression

params <- fitdistr(x = log(gene_expr), "normal")
params[[1]]
nm_mean <- params[[1]][1]
nm_sd <- params[[1]][2]


ggp <- ggplot(data = data.frame(gene_expr = gene_expr), aes(x = log(gene_expr))) + 
  geom_density() +
  geom_line(data = data.frame(x = seq(min(log(gene_expr)), max(log(gene_expr)), by = 0.01), y = dnorm(seq(min(log(gene_expr)), max(log(gene_expr)), by = 0.01), mean = nm_mean, sd = nm_sd)), aes(x = x, y = y), colour = "red") +
  ggtitle(paste0("Normal fit with mean = ", round(nm_mean, 2), " and sd = ", round(nm_sd, 2)))

pdf(paste0(out_dir, "nm_",data_name,"_", count_method, "_hist_normal.pdf"))
print(ggp)
dev.off()


write.table(round(params[[1]], 2), file = paste0(out_dir, "nm_",data_name,"_", count_method, "_normal.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)





##############################################################################
### proportions exact as in real data
##############################################################################

### calculate the distribution of proportions from one sample


metadata <- read.table(paste0(main_data_dir, "/3_metadata/metadata.xls"), stringsAsFactors = FALSE, sep="\t", header=TRUE) 

metadata_org <- metadata


count_dir <- paste0(main_data_dir, "/2_counts/", count_method, "/")


### load counts
counts_list <- lapply(1:length(metadata_org$sampleName), function(i){
  # i = 1
  cts <- read.table(paste0(count_dir, metadata_org$sampleName[i], ".counts"), header = FALSE, as.is = TRUE)
  colnames(cts) <- c("group_id", metadata_org$sampleName[i])  
  return(cts)
})

counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)
counts <- counts[!grepl(pattern = "_", counts$group_id),]

group_split <- strsplit2(counts[,1], ":")
counts <- counts[, -1]



d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = colnames(counts), group = rep("C1", ncol(counts)))

d <- d[, use_sample]
lib_size <- sum(counts[ , use_sample])*1e-6 ### for cpm calculations

d <- dmFilter(d, min_samps_gene_expr = 1, min_samps_feature_prop = 0, min_feature_prop = 0, min_gene_expr = 100/lib_size)

plotData(d, out_dir = paste0(out_dir, "prop_", data_name, "_", count_method, "_orig_"))





calculate_proportions <- function(d, name_approach = "", out_dir, data_name, count_method, nr_features = seq(2, 15, 1)){
  
  
  cts <- d@counts[, 1]
  
  
  prop_list <- lapply(1:length(cts), function(g){
    # g = 1
    
    if(sum(cts[[g]]) == 0)
      return(NULL)
    
    pi <- sort(cts[[g]]/sum(cts[[g]]), decreasing = TRUE)
    
    out <- data.frame(gene_id = paste0("g", g), feature_id = paste0("f", 1:length(pi)), proportions = pi)
    
    return(out)
    
  })
  
  
  prop <- rbind.fill(prop_list)
  
  write.table(prop, file = paste0(out_dir, "prop_",data_name,"_", count_method, name_approach, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
  
  ##############################################################################
  ### proportions per number of features 
  ##############################################################################
  
  ### Uniform 
  
  write.table(rep(1, 3)/3, file = paste0(out_dir, "prop_q3_uniform.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  write.table(rep(1, 10)/10, file = paste0(out_dir, "prop_q10_uniform.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  
  
  ### Decay from real data 
  
  max_features <- max(width(cts))
  max_features
  max_features_plot <- min(15, max_features)
  nr_features <- nr_features[nr_features <= max_features]
  
  
  prop_list <- lapply(1:length(cts), function(g){
    # g = 1
    if(sum(cts[[g]]) == 0)
      return(NULL)
    
    pi <- c(sort(cts[[g]]/sum(cts[[g]]), decreasing = TRUE), rep(NA, max_features + 1 - nrow(cts[[g]])))
    pi[max_features + 1] <- nrow(cts[[g]])
    
    return(pi)
  })
  
  
  prop <- do.call(rbind, prop_list)
  prop <- data.frame(prop)
  colnames(prop) <- c(paste0("F", 1:max_features), "Nr_features")
  
  
  propm <- melt(prop, id.vars = "Nr_features", variable.name = "Features", value.name = "Proportions") 
  propm <- propm[complete.cases(propm), ]
  propm$Nr_features <- factor(propm$Nr_features)
  
  
  pdf(paste0(out_dir, "prop_",data_name,"_", count_method, name_approach, "_boxplots_overall.pdf"), 15, 5)
  ggp <- ggplot(propm, aes(x = Features, y = Proportions)) +
    geom_boxplot(outlier.size = 0.5, outlier.colour = "brown") +
    xlab("Sorted features")
  print(ggp)
  dev.off()
  
  
  propm_sub <- propm[as.numeric(as.character(propm$Nr_features)) <= max_features_plot, ]
  propm_sub$Nr_features <- factor(propm_sub$Nr_features)
  propm_sub$Features <- factor(propm_sub$Features)
  levels(propm_sub$Nr_features) <- paste0(levels(propm_sub$Nr_features), " (", as.numeric(table(propm_sub$Nr_features)), ")")
  levels(propm_sub$Features)
  
  
  pdf(paste0(out_dir, "prop_",data_name,"_", count_method, name_approach, "_boxplots.pdf"), width = 15)
  ggp <- ggplot(propm_sub, aes(x = Features, y = Proportions, fill = Nr_features)) +
    geom_boxplot() +
    xlab("Sorted features") +
    scale_fill_discrete(name = "Total number \nof features") +
    coord_cartesian(ylim = c(-0.1, 1.1)) 
  print(ggp)
  dev.off()
  
  
  
  ### generate proportions from proportions per Nr_features
  
  gen_prop_list <- lapply(nr_features, function(i){
    # i = 19
    # print(i)
    prop_tmp <- prop[prop$Nr_features == i, 1:i]
    
    if(nrow(prop_tmp) == 0)
      return(NULL)
    
    prop_dir <- colMedians(as.matrix(prop_tmp))
    prop_dir <- sort(prop_dir/sum(prop_dir), decreasing = TRUE)
    # print(sum(colMedians(as.matrix(prop_tmp))))
    
    write.table(prop_dir, file = paste0(out_dir, "prop_q", i, "_",data_name,"_", count_method, name_approach, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    out <- c(prop_dir, rep(NA, max_features + 1 - i))
    out[max_features + 1] <- i
    
    return(out)
    
  })
  
  
  gen_prop <- do.call(rbind, gen_prop_list)
  gen_prop <- data.frame(gen_prop)
  colnames(gen_prop) <- c(paste0("F", 1:max_features), "Nr_features")
  
  gen_propm <- melt(gen_prop, id.vars = "Nr_features", variable.name = "Features", value.name = "Proportions") 
  gen_propm <- gen_propm[complete.cases(gen_propm), ]
  gen_propm$Nr_features <- factor(gen_propm$Nr_features)
  
  gen_propm_sub <- gen_propm[as.numeric(as.character(gen_propm$Nr_features)) <= max_features_plot, ]
  gen_propm_sub$Nr_features <- factor(gen_propm_sub$Nr_features)
  gen_propm_sub$Features <- factor(gen_propm_sub$Features)
  levels(gen_propm_sub$Nr_features)
  levels(gen_propm_sub$Features)
  
  
  pdf(paste0(out_dir, "prop_",data_name,"_", count_method, name_approach, "_parameters.pdf"), width = 15)
  ggp <- ggplot(gen_propm_sub, aes(x = Features, y = Proportions, group = Nr_features, colour = Nr_features)) +
    geom_line() +
    xlab("Sorted features") +
    scale_colour_discrete(name = "Total number \nof features") +
    coord_cartesian(ylim = c(-0.1, 1.1))
  
  print(ggp)
  dev.off()
  
  
  
  ### generate proportions from overall proportions
  
  gen_prop_list <- lapply(nr_features, function(i){
    # i = 19
    # print(i)
    
    prop_dir <- colMedians(as.matrix(prop[, 1:max_features]), na.rm = TRUE)[1:i]
    prop_dir <- sort(prop_dir/sum(prop_dir), decreasing = TRUE)
    # print(sum(colMedians(as.matrix(prop_tmp))))
    
    write.table(prop_dir, file = paste0(out_dir, "prop_q", i, "_",data_name,"_", count_method, name_approach, "_overall.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    out <- c(prop_dir, rep(NA, max_features + 1 - i))
    out[max_features + 1] <- i
    
    return(out)
    
  })
  
  
  gen_prop <- do.call(rbind, gen_prop_list)
  gen_prop <- data.frame(gen_prop)
  colnames(gen_prop) <- c(paste0("F", 1:max_features), "Nr_features")
  
  gen_propm <- melt(gen_prop, id.vars = "Nr_features", variable.name = "Features", value.name = "Proportions") 
  gen_propm <- gen_propm[complete.cases(gen_propm), ]
  gen_propm$Nr_features <- factor(gen_propm$Nr_features)
  
  gen_propm_sub <- gen_propm[as.numeric(as.character(gen_propm$Nr_features)) <= max_features_plot, ]
  gen_propm_sub$Nr_features <- factor(gen_propm_sub$Nr_features)
  gen_propm_sub$Features <- factor(gen_propm_sub$Features)
  levels(gen_propm_sub$Nr_features)
  levels(gen_propm_sub$Features)
  
  
  
  gen_prop_overall <- colMedians(as.matrix(prop[, 1:max_features]), na.rm = TRUE)
  gen_prop_overall <- data.frame(Nr_features = 2, Features = paste0("F", 1:max_features), Proportions = gen_prop_overall/sum(gen_prop_overall))
  
  gen_prop_overall_sub <- gen_prop_overall[1:max_features_plot, ]
  gen_prop_overall_sub$Nr_features <- factor(gen_prop_overall_sub$Nr_features, levels = levels(gen_propm_sub$Nr_features))
  gen_prop_overall_sub$Features <- factor(gen_prop_overall_sub$Features, levels = levels(gen_propm_sub$Features)
  )
  levels(gen_prop_overall_sub$Nr_features)
  levels(gen_prop_overall_sub$Features)
  
  
  
  pdf(paste0(out_dir, "prop_",data_name,"_", count_method, name_approach, "_parameters_overall.pdf"), width = 15)
  ggp <- ggplot(gen_propm_sub, aes(x = Features, y = Proportions, group = Nr_features, colour = Nr_features)) +
    geom_line() +
    xlab("Sorted features") +
    scale_colour_discrete(name = "Total number \nof features") +
    geom_line(data = gen_prop_overall_sub, aes(x = Features, y = Proportions), colour = "black") +
    coord_cartesian(ylim = c(-0.1, 1.1))
  
  print(ggp)
  dev.off()
  
  
  
}




### Approach with using raw counts for one sample and filtering on it with dmFilter - keep features with min_feature_prop = 0.01


d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = colnames(counts), group = rep("C1", ncol(counts)))

d <- d[, use_sample]

lib_size <- sum(counts[ , use_sample])*1e-6 ### for cpm calculations

d <- dmFilter(d, min_samps_gene_expr = 1, min_samps_feature_prop = 1, min_feature_prop = 0.01, min_gene_expr = 100/lib_size)


plotData(d, out_dir = paste0(out_dir, "prop_", data_name, "_", count_method, "_"))


calculate_proportions(d, name_approach = "", out_dir, data_name, count_method, nr_features = seq(2, 15, 1))




### Approach with filtering features based on their expression (in CPM)

lib_size <- sum(counts[, use_sample])*1e-6
feature_cutoff <- round(lib_size/2)
feature_cutoff

keep_counts <- counts[, use_sample] > feature_cutoff

counts_filt <- counts[keep_counts, , drop = FALSE] 
group_split_filt <- group_split[keep_counts, , drop = FALSE]

lib_size <- sum(counts_filt)*1e-6


d <- dmDSdata(counts = counts_filt, gene_id = group_split_filt[, 1], feature_id = group_split_filt[, 2], sample_id = colnames(counts_filt), group = rep("C1", ncol(counts_filt)))

d <- d[, use_sample]

d <- dmFilter(d, min_samps_gene_expr = 1, min_samps_feature_prop = 1, min_feature_prop = 0, min_gene_expr = 100/lib_size)


plotData(d, out_dir = paste0(out_dir, "prop_", data_name, "_", count_method, "_fcutoff_"))


calculate_proportions(d, name_approach = paste0("_fcutoff"), out_dir, data_name, count_method, nr_features = c(seq(2, 15, 1), 20, 25, 30, 35, 40, 45, 50))






















