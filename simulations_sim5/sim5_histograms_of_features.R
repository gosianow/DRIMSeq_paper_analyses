######################################################
## ----- sim5_histograms_of_features
## <<sim5_histograms_of_features.R>>

# BioC 3.2
# Created 19 Nov 2015 
# Modified 14 Dec 2015

##############################################################################

library(iCOBRA)
library(Hmisc)
library(DEXSeq)
library(DRIMSeq)
library(limma)
library(plyr)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_sim5'
simulation='hsapiens_node_nonull'
count_method=c('htseq','kallisto','htseq_prefiltered15')[2]
filter_method='filter1'

##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(simulation)
print(count_method)
print(filter_method)


##############################################################################

setwd(paste0(rwd, "/", simulation))

method_out <- "drimseq_0_3_3"
comparison_out <- "drimseq_0_3_3_comparison"
out_dir <- paste0(comparison_out, "/", filter_method, "/", count_method, "_")

dir.create(dirname(out_dir), recursive = TRUE, showWarnings = FALSE)


#######################################################
# DEXSeq results
#######################################################


if(count_method == "htseq")
  results_dir <- "4_results/dexseq_htseq_nomerge"
if(count_method == "kallisto")
  results_dir <- "4_results/dexseq_kallisto"
if(count_method == "htseqprefiltered15")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast15"
if(count_method == "htseqprefiltered5")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast5"
if(count_method == "kallistofiltered5")
  results_dir <- "4_results/dexseq_kallisto_txfilt_5"
if(count_method == "kallistoprefiltered5")
  results_dir <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_kallisto_kallistoest_atleast5"

results_dir


### DEXSeq exon results
load(paste0(results_dir, ".Rdata"))
# res

rt <- as.data.frame(res[, 1:7])

table(is.na(rt$dispersion))
table(is.na(rt$pvalue))
table(is.na(rt$padj))


rt <- rt[, c("groupID", "featureID", "dispersion", "pvalue")]

table(table(rt$groupID), useNA = "always")


keep <- complete.cases(rt)
table(keep)

rt <- rt[keep, ]

table(table(rt$groupID), useNA = "always")


### Gene level results

# rtg <- read.table(paste0(results_dir, ".txt"), header = TRUE, as.is = TRUE)
# colnames(rtg) <- c("gene_id", "dexseq")
# dim(rtg)
# rtg <- rtg[complete.cases(rtg), ]
# dim(rtg)
# 
# all(names(table(rt$groupID)) %in% rtg$gene_id)



#################################################################
# Combined plot 
#################################################################

results_dir <- paste0(method_out, "/", count_method, "/", filter_method, "/")

### All

load(paste0(results_dir, "d.Rdata"))

d_filt <- dmFilter(d, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0, max_features = Inf)

df_all <- data.frame(nr_features = width(d_filt@counts), Genes = paste0("all (", length(width(d_filt@counts)), ")"))


### DRIMSeq

load(paste0(results_dir, "drimseq_genewise_grid_none_d.Rdata"))

df_drimseq <- data.frame(nr_features = width(d@counts), Genes = paste0("drimseq (", length(width(d@counts)), ")"))


### DEXSeq
df_dexseq <- data.frame(nr_features = as.numeric(table(rt$groupID)), Genes = paste0("dexseq (", length(as.numeric(table(rt$groupID))), ")"))

df <- rbind.fill(df_all, df_dexseq, df_drimseq)

levels(df$Genes)


ggp <- ggplot(df, aes(x = nr_features, y = ..count.., colour = Genes, linetype = Genes)) +
  geom_freqpoly(binwidth = 1, size = 1.3, alpha = 0.6) +
  theme_bw() +
  xlab("Number of features per gene") +
  ylab("Frequency") +
  # coord_cartesian(xlim = c(0, max(df_drimseq$nr_features))) +
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 16, face = "bold"), legend.position = "bottom", strip.text = element_text(size = 14)) +
  scale_linetype_manual(values = c(2, 1, 1))


pdf(paste0(out_dir, "hist_features.pdf"), 7, 5)
print(ggp)
dev.off()




# ggp <- ggplot(df, aes(x = nr_features, y = ..count.., colour = Genes)) +
#   geom_line(stat="density", size = 1.2, alpha = 1, adjust = 2, linetype = 1) +
#   theme_bw() +
#   xlab("Log2 of number of features per gene") +
#   ylab("Frequency") +
#   theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 16, face = "bold"), legend.position = "bottom", strip.text = element_text(size = 14)) 
# 
# 
# pdf(paste0(out_dir, "hist_features_density.pdf"), 7, 5)
# print(ggp)
# dev.off()












