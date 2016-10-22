######################################################
## <<sim5_dexseq_run.R>>

# BioC 2.14

# Created 15 Nov 2015 
# Updated 18 Oct 2016

##############################################################################

library("BiocParallel")
library("DEXSeq")

##############################################################################
# Test arguments
##############################################################################

rwd='/home/gosia/multinomial_project/simulations_sim5/hsapiens_withde_nonull'
workers=4
count_method=c('htseq','kallisto','htseqprefiltered5','kallistofiltered5','kallistoprefiltered5')[5]

##############################################################################
# Read in the arguments
##############################################################################

rm(list = ls())

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

cat(paste0(args, collapse = "\n"), fill = TRUE)

##############################################################################

setwd(rwd)

if(workers > 1){
  BPPARAM <- MulticoreParam(workers = workers)
}else{
  BPPARAM <- SerialParam()
}


if(count_method == "htseq"){
  count_dir <- "2_counts/dexseq_nomerge/dexseq"
  method_name <- "dexseq_htseq_nomerge"
  out_basename <- "4_results/dexseq_htseq_nomerge"
}
if(count_method == "kallisto"){
  count_dir <- "2_counts/kallisto/kallisto"
  method_name <- "dexseq_kallisto"
  out_basename <- "4_results/dexseq_kallisto"
}
if(count_method == "htseqprefiltered5"){
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast5/dexseq"
  method_name <- "dexseq_htseq_nomerge_kallistoest_atleast5"
  out_basename <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_htseq_nomerge_kallistoest_atleast5"
}
if(count_method == "kallistofiltered5"){
  count_dir <- "2_counts/kallisto_txfilt_5/kallisto"
  method_name <- "dexseq_kallisto_txfilt_5"
  out_basename <- "4_results/dexseq_kallisto_txfilt_5"
}
if(count_method == "kallistoprefiltered5"){
  count_dir <- "2_counts/INCOMPLETE_KALLISTOEST/kallisto_kallistoest_atleast5/kallisto"
  method_name <- "dexseq_kallisto_kallistoest_atleast5"
  out_basename <- "4_results/INCOMPLETE_KALLISTOEST/dexseq_kallisto_kallistoest_atleast5"
}

count_dir
method_name
out_basename

dir.create("4_results/INCOMPLETE_KALLISTOEST", recursive = TRUE)


###############################################################################
# metadata
###############################################################################


metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata


###############################################################################
# Run DEXSeq
###############################################################################

conditions <- metadata$group

## Generate sample table
sample_table <- data.frame(sample = 1:length(conditions), condition = conditions)

## Generate DEXSeqDataSet

count_files <- paste0(count_dir, 1:length(conditions), ".txt")

dxd <- DEXSeqDataSetFromHTSeq(count_files, sampleData = sample_table, design = ~ sample + exon + condition:exon, flattenedfile = NULL)

## Estimate size factors and dispersions
dxd <- estimateSizeFactors(dxd)
message("Estimating dispersions...")
dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)

## Run the test and get results for exon bins and genes
message("Running test...")
dxd <- testForDEU(dxd, BPPARAM = BPPARAM)

message("Summarizing results on gene level...")
res <- DEXSeqResults(dxd)
pgq <- perGeneQValue(res, p = "pvalue")

## Save results
save(dxd, res, pgq, file = paste0(out_basename, ".Rdata"))


tmp <- cbind(gene = names(pgq), "adjP" = pgq)
colnames(tmp)[colnames(tmp) == "adjP"] <- paste0(method_name, ":adjP")

write.table(tmp, file = paste0(out_basename, ".txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")



sessionInfo()















##################################
### sim5_dexseq_run.R done!
##################################