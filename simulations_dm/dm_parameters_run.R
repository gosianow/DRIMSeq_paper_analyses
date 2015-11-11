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


##############################################################################
# Read in the arguments
##############################################################################

# rwd='/home/gosia/multinomial_project/simulations_dm/drimseq_0_3_1/'
# count_method=c('htseq','kallisto')[2]


## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(count_method)



##############################################################################
# Parameters based on Kim data
##############################################################################


setwd(rwd)

out_dir <- "dm_parameters/"
dir.create(out_dir, recursive = T, showWarnings = FALSE)

main_data_dir <- "/home/Shared/data/seq/kim_adenocarcinoma/"
method_out <- "drimseq_0_3_1"
model <- "model_full"

data_dir <- paste0(main_data_dir, method_out, "/",  model, "/", count_method, "/")

disp <- "genewise"
disp_mode <- "grid"
disp_moderation <- "none"

data_name <- paste0(data_dir, "/drimseq_", disp, "_", disp_mode, "_", disp_moderation, "_")

load(paste0(data_name, "d.Rdata"))


#######################################
### common dispersion
#######################################

common_disp <- common_dispersion(d)

write.table(round(common_disp), file = paste0(out_dir, "disp_common_kim_", count_method, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#######################################
### genewise dispersion
#######################################

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



#######################################
### dispersion for gene expression NB()
#######################################


gene_counts <- lapply(1:length(d@counts), function(g){ colSums(d@counts[[g]]) })
gene_counts <- do.call(rbind, gene_counts)
rownames(gene_counts) <- names(d)


dge <- DGEList(gene_counts, group = samples(d)$group)

dge <- calcNormFactors(dge)

pdf(paste0(out_dir, "gene_expr_kim_", count_method, "_MDS.pdf"))
plotMDS(dge, col = c("darkgreen", "blue")[samples(d)$group], xlim = c(-2, 2), ylim = c(-2, 2))
dev.off()


dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

pdf(paste0(out_dir, "gene_expr_kim_", count_method, "_BCV.pdf"))
plotBCV(dge)
dev.off()

dge$common.dispersion


write.table(round(dge$common.dispersion, 2), file = paste0(out_dir, "gene_expr_disp_common_kim_", count_method, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


bcv <- sqrt(dge$common.dispersion)
bcv



pdf(paste0(out_dir, "gene_expr_disp_tagwise_kim_", count_method, "_hist.pdf"))
ggp <- ggplot(data = data.frame(tagwise_dispersion = dge$tagwise.dispersion), aes(x = tagwise_dispersion)) + 
  geom_density() 
print(ggp)
dev.off()




#######################################
### proportions
#######################################

write.table(rep(1, 3)/3, file = paste0(out_dir, "prop_q3_uniform.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(rep(1, 10)/10, file = paste0(out_dir, "prop_q10_uniform.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



### calculate the distribution of proportions from one sample

# ### Approach with filtering features based on their CPM
# counts <- read.table(paste0(main_data_dir, "2_counts/", count_method, "/", count_method, "_counts.txt"), header = TRUE, as.is = TRUE)
# counts <- counts[!grepl(pattern = "_", counts$group_id),]
# 
# lib_size <- sum(counts[, "GSM927308"])*1e-6
# 
# counts <- counts[counts[, "GSM927308"] > 10, , drop = FALSE] # round(lib_size/2)
# 
# lib_size <- sum(counts[, "GSM927308"])*1e-6
# 
# group_split <- strsplit2(counts[,1], ":")
# counts <- counts[, -1]
# 
# d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = colnames(counts), group = rep("C1", ncol(counts)))
# 
# d <- d[, "GSM927308"]
# 
# d <- dmFilter(d, min_samps_gene_expr = 1, min_samps_feature_prop = 1, min_feature_prop = 0, min_gene_expr = 100/lib_size)
# # d <- dmFilter(d, min_samps_gene_expr = 1, min_samps_feature_prop = 1, min_feature_prop = 0.01, min_gene_expr = 100/lib_size)
# 
# plotData(d, out_dir = paste0(out_dir, "prop_kim_", count_method, "_gene_expression_"))
# 
# cts <- d@counts[, 1]
# 
# gene_expr <- sapply(1:length(cts), function(g){ sum(cts[[g]]) })
# 
# pdf(paste0(out_dir, "prop_kim_", count_method, "_gene_expression.pdf"))
# ggp <- ggplot(data.frame(gene_expr = gene_expr), aes(x = log10(gene_expr))) +
#   geom_density()
# print(ggp)
# dev.off()




### Approach with using raw counts for one sample and filtering on it with dmFilter
counts <- read.table(paste0(main_data_dir, "2_counts/", count_method, "/", count_method, "_counts.txt"), header = TRUE, as.is = TRUE)
counts <- counts[!grepl(pattern = "_", counts$group_id),]
group_split <- strsplit2(counts[,1], ":")
counts <- counts[, -1]

d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = colnames(counts), group = rep("C1", ncol(counts)))

d <- d[, "GSM927308"]

lib_size <- sum(counts[, "GSM927308"])*1e-6

d <- dmFilter(d, min_samps_gene_expr = 1, min_samps_feature_prop = 1, min_feature_prop = 0.01, min_gene_expr = 100/lib_size)
plotData(d, out_dir = paste0(out_dir, "prop_kim_", count_method, "_gene_expression_"))

cts <- d@counts[, 1]
gene_expr <- sapply(1:length(cts), function(g){ sum(cts[[g]]) })

pdf(paste0(out_dir, "prop_kim_", count_method, "_gene_expression.pdf"))
ggp <- ggplot(data.frame(gene_expr = gene_expr), aes(x = log10(gene_expr))) +
  geom_density()
print(ggp)
dev.off()





### Approach with using d from the analysis
# d <- d[, "GSM927308"]
# plotData(d, out_dir = paste0(out_dir, "prop_kim_", count_method, "_gene_expression_"))
# cts <- d@counts[, 1]
# gene_expr <- sapply(1:length(cts), function(g){ sum(cts[[g]]) })
# pdf(paste0(out_dir, "prop_kim_", count_method, "_gene_expression.pdf"))
# ggp <- ggplot(data.frame(gene_expr = gene_expr), aes(x = log10(gene_expr))) +
#   geom_density()
# print(ggp)
# dev.off()
# cts <- cts[gene_expr > 500, ]



max_features <- max(width(cts))
max_features


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


pdf(paste0(out_dir, "prop_kim_", count_method, "_boxplots_overall.pdf"))
ggp <- ggplot(propm, aes(x = Features, y = Proportions)) +
  geom_boxplot() +
  xlab("Sorted features")
print(ggp)
dev.off()


propm_sub <- propm[as.numeric(as.character(propm$Nr_features)) <= 15, ]
propm_sub$Nr_features <- factor(propm_sub$Nr_features)
propm_sub$Features <- factor(propm_sub$Features)
levels(propm_sub$Nr_features) <- paste0(levels(propm_sub$Nr_features), " (", as.numeric(table(propm_sub$Nr_features)), ")")
levels(propm_sub$Features)


pdf(paste0(out_dir, "prop_kim_", count_method, "_boxplots.pdf"), width = 15)
ggp <- ggplot(propm_sub, aes(x = Features, y = Proportions, fill = Nr_features)) +
  geom_boxplot() +
  xlab("Sorted features") +
  scale_fill_discrete(name = "Total number \nof features") +
  coord_cartesian(ylim = c(-0.1, 1.1)) 
print(ggp)
dev.off()



### generate proportions from proportions per Nr_features

nr_features <- seq(2, 15, 1)

gen_prop_list <- lapply(nr_features, function(i){
  # i = 19
  # print(i)
  prop_tmp <- prop[prop$Nr_features == i, 1:i]
  
  prop_dir <- colMedians(as.matrix(prop_tmp))
  prop_dir <- prop_dir/sum(prop_dir)
  # print(sum(colMedians(as.matrix(prop_tmp))))

  write.table(prop_dir, file = paste0(out_dir, "prop_q", i, "_kim_", count_method, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
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

gen_propm_sub <- gen_propm[as.numeric(as.character(gen_propm$Nr_features)) <= 15, ]
gen_propm_sub$Nr_features <- factor(gen_propm_sub$Nr_features)
gen_propm_sub$Features <- factor(gen_propm_sub$Features)
levels(gen_propm_sub$Nr_features)
levels(gen_propm_sub$Features)


pdf(paste0(out_dir, "prop_kim_", count_method, "_parameters.pdf"), width = 15)
ggp <- ggplot(gen_propm_sub, aes(x = Features, y = Proportions, group = Nr_features, colour = Nr_features)) +
  geom_line() +
  xlab("Sorted features") +
  scale_colour_discrete(name = "Total number \nof features") +
  coord_cartesian(ylim = c(-0.1, 1.1))

print(ggp)
dev.off()







### generate proportions from overall proportions

nr_features <- seq(2, 15, 1)

gen_prop_list <- lapply(nr_features, function(i){
  # i = 19
  # print(i)

  prop_dir <- colMedians(as.matrix(prop[, 1:max_features]), na.rm = TRUE)[1:i]
  prop_dir <- prop_dir/sum(prop_dir)
  # print(sum(colMedians(as.matrix(prop_tmp))))
  
  write.table(prop_dir, file = paste0(out_dir, "prop_q", i, "_kim_", count_method, "_overall.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
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

gen_propm_sub <- gen_propm[as.numeric(as.character(gen_propm$Nr_features)) <= 15, ]
gen_propm_sub$Nr_features <- factor(gen_propm_sub$Nr_features)
gen_propm_sub$Features <- factor(gen_propm_sub$Features)
levels(gen_propm_sub$Nr_features)
levels(gen_propm_sub$Features)



gen_prop_overall <- colMedians(as.matrix(prop[, 1:max_features]), na.rm = TRUE)
gen_prop_overall <- data.frame(Nr_features = 2, Features = paste0("F", 1:max_features), Proportions = gen_prop_overall/sum(gen_prop_overall))

gen_prop_overall_sub <- gen_prop_overall[1:15, ]
gen_prop_overall_sub$Nr_features <- factor(gen_prop_overall_sub$Nr_features, levels = levels(gen_propm_sub$Nr_features))
gen_prop_overall_sub$Features <- factor(gen_prop_overall_sub$Features, levels = levels(gen_propm_sub$Features)
)
levels(gen_prop_overall_sub$Nr_features)
levels(gen_prop_overall_sub$Features)



pdf(paste0(out_dir, "prop_kim_", count_method, "_parameters_overall.pdf"), width = 15)
ggp <- ggplot(gen_propm_sub, aes(x = Features, y = Proportions, group = Nr_features, colour = Nr_features)) +
  geom_line() +
  xlab("Sorted features") +
  scale_colour_discrete(name = "Total number \nof features") +
  geom_line(data = gen_prop_overall_sub, aes(x = Features, y = Proportions), colour = "black") +
  coord_cartesian(ylim = c(-0.1, 1.1))

print(ggp)
dev.off()





























