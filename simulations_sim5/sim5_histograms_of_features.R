######################################################
## ----- sim5_histograms_of_features
## <<sim5_histograms_of_features.R>>

# BioC 3.1
# Created 19 Nov 2015 

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



##############################################################################


setwd(paste0(rwd, "/", simulation))

########################################################
# create metadata file
########################################################

metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata

#######################################################
# plot histograms of features for DEXSeq
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




tt <- as.numeric(table(rt$groupID))

df <- data.frame(tt = tt)
binwidth <- ceiling(max(df$tt)/50)

ggp <- ggplot(df, aes_string(x = "tt")) +
  theme_bw() +
  xlab("Number of features per gene") +
  ylab("Frequency") +
  geom_histogram(fill = "olivedrab", binwidth = binwidth) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold")) +
  coord_cartesian(xlim = c(0, max(tt) + 2)) +
  geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(length(tt), " genes   \n ", sum(tt) , " features   ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, size = 6)


pdf(paste0(results_dir, "_orig_data_all_hist_features.pdf"))
print(ggp)
dev.off()





tt <- tt[tt > 1]

df <- data.frame(tt = tt)
binwidth <- ceiling(max(df$tt)/50)

ggp <- ggplot(df, aes_string(x = "tt")) +
  theme_bw() +
  xlab("Number of features per gene") +
  ylab("Frequency") +
  geom_histogram(fill = "olivedrab", binwidth = binwidth) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold")) +
  coord_cartesian(xlim = c(0, max(tt) + 2)) +
  geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(length(tt), " genes   \n ", sum(tt) , " features   ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, size = 6)


pdf(paste0(results_dir, "_orig_data_no1feature_hist_features.pdf"))
print(ggp)
dev.off()




keep <- complete.cases(rt)
table(keep)
rt <- rt[keep, ]


table(table(rt$groupID), useNA = "always")

table(as.numeric(table(rt$groupID)) > 1)
table(as.numeric(table(rt$groupID)) > 0)



tt <- as.numeric(table(rt$groupID))

df <- data.frame(tt = tt)
binwidth <- ceiling(max(df$tt)/50)

ggp <- ggplot(df, aes_string(x = "tt")) +
  theme_bw() +
  xlab("Number of features per gene") +
  ylab("Frequency") +
  geom_histogram(fill = "olivedrab", binwidth = binwidth) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold")) +
  coord_cartesian(xlim = c(0, max(tt) + 2)) +
  geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(length(tt), " genes   \n ", sum(tt) , " features   ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, size = 6)


pdf(paste0(results_dir, "_hist_features.pdf"))
print(ggp)
dev.off()



### Gene level results

# rtg <- read.table(paste0(results_dir, ".txt"), header = TRUE, as.is = TRUE)
# colnames(rtg) <- c("gene_id", "dexseq")
# dim(rtg)
# rtg <- rtg[complete.cases(rtg), ]
# dim(rtg)

# all(names(table(rt$groupID)) %in% rtg$gene_id)


#######################################################
# plot histograms of features for DRIMSeq
#######################################################


if(count_method == "htseq")
  count_dir <- "2_counts/dexseq_nomerge/dexseq"
if(count_method == "kallisto")
  count_dir <- "2_counts/kallisto/kallisto"
if(count_method == "htseqprefiltered15")
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast15/dexseq"
if(count_method == "htseqprefiltered5")
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast5/dexseq"
if(count_method == "kallistofiltered5")
  count_dir <- "2_counts/kallisto_txfilt_5/kallisto"
if(count_method == "kallistoprefiltered5")
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/kallisto_kallistoest_atleast5/kallisto"

count_dir

### load counts
counts_list <- lapply(1:6, function(i){
  # i = 1
  cts <- read.table(paste0(count_dir, i, ".txt"), header = FALSE, as.is = TRUE)
  colnames(cts) <- c("group_id", paste0("sample_", i))  
  return(cts)
})

counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)

counts <- counts[!grepl(pattern = "_", counts$group_id),]
group_split <- strsplit2(counts[,1], ":")
counts <- counts[, -1]




method_out <- "drimseq_0_3_1"
out_dir <- paste0(method_out, "/", count_method, "/")


d_org <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = metadata$sample_id, group = metadata$group)

plotData(d_org, out_dir = paste0(out_dir, "orig_data_all_"))



d_filt <- dmFilter(d_org, min_samps_gene_expr = 0, min_samps_feature_prop = 0, min_gene_expr = 0, min_feature_prop = 0)

plotData(d_filt, out_dir = paste0(out_dir, "orig_data_no1f_no0g_"))


#################################################################
# Combined plot 
#################################################################

comparison_out <- "drimseq_0_3_1_comparison"
out_dir <- paste0(comparison_out, "/", count_method, "_")

### DRIMSeq
df_all <- data.frame(nr_features = width(d_filt@counts), Genes = paste0("all (", length(width(d_filt@counts)), ")"))

load(paste0(method_out, "/", count_method, "/d.Rdata"))

df_drimseq <- data.frame(nr_features = width(d@counts), Genes = paste0("drimseq (", length(width(d@counts)), ")"))


### DEXSeq
df_dexseq <- data.frame(nr_features = as.numeric(table(rt$groupID)), Genes = paste0("dexseq (", length(as.numeric(table(rt$groupID))), ")"))

df <- rbind.fill(df_all, df_dexseq, df_drimseq)





ggp <- ggplot(df, aes(x = nr_features, y = ..count.., colour = Genes)) +
  theme_bw() +
  xlab("Number of features per gene") +
  ylab("Frequency") +
  geom_freqpoly(binwidth = 1, size = 1.5, alpha = 0.8) +
  # coord_cartesian(xlim = c(0, max(df_drimseq$nr_features))) +
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 16, face = "bold"), legend.position = "bottom", strip.text = element_text(size = 14)) 


pdf(paste0(out_dir, "hist_features.pdf"), 7, 5)
print(ggp)
dev.off()




ggp <- ggplot(df, aes(x = nr_features, y = ..count.., colour = Genes)) +
  theme_bw() +
  xlab("Log2 of number of features per gene") +
  ylab("Frequency") +
  geom_line(stat="density", size = 1.2, alpha = 1, adjust = 2, linetype = 1) +
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), axis.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 16, face = "bold"), legend.position = "bottom", strip.text = element_text(size = 14)) 


pdf(paste0(out_dir, "hist_features_density.pdf"), 7, 5)
print(ggp)
dev.off()












