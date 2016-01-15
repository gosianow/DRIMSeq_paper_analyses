######################################################
## ----- kim_kallisto_filter_gtf
## <<kim_kallisto_filter_gtf.R>>

# BioC 3.1
# Crated 24 Nov 2015
# Get kallisto counts for reduced annotation

#######################################################

Sys.time()

#######################################################

library(DRIMSeq)
library(limma)
library(rtracklayer)
library(tools)

##############################################################################
# Test arguments
##############################################################################


# rwd='/home/Shared/data/seq/kim_adenocarcinoma'
# gtf_path='/home/Shared/data/annotation/Human/Ensembl_GRCh37.71/gtf/Homo_sapiens.GRCh37.71.gtf'



##############################################################################
# Read in the arguments
##############################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)



##############################################################################


setwd(rwd)

gtf_new_path <- paste0(dirname(gtf_path) ,"/", basename(file_path_sans_ext(gtf_path)) ,"_kallistoest_atleast5.gtf")



##########################################################################
# load metadata
##########################################################################


metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors = FALSE, sep="\t", header=TRUE) 


#########################################################################################
# Use TPMs fractions to filter out the low expressed transcripts from 
# - gtf annotation
# - kallisto expected counts
#########################################################################################


gtf <- import(gtf_path)

gene_trans <- unique(mcols(gtf)[, c("gene_id", "transcript_id")])
rownames(gene_trans) <- gene_trans$transcript_id


### load TPMs

samples <- metadata[, "sampleName"]

tpm_list <- lapply(1:length(samples), function(i){
  # i = 1
  print(i)
  
  abundance <- read.table(paste0(rwd, "2_counts/kallisto/", samples[i], "/", "abundance.txt"), header = TRUE, sep = "\t", as.is = TRUE)
  
  tpm <- data.frame(paste0(gene_trans[abundance[, "target_id"], "gene_id"], ":", abundance[, "target_id"]), tpm = abundance[, "tpm"], stringsAsFactors = FALSE)
  colnames(tpm) <- c("group_id", samples[i])
  
  return(tpm)
  
})



tpm <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), tpm_list)

group_split <- strsplit2(tpm[, "group_id"], ":")

d_org <- dmDSdata(counts = tpm[, -1], gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = colnames(tpm[, -1]), group = rep("C1", ncol(tpm[, -1])))

d_filt <- dmFilter(d_org, min_samps_gene_expr = 0, min_samps_feature_prop = 1, min_gene_expr = 0, min_feature_prop = 0.05)

trans_keep <- counts(d_filt)$feature_id



######################## - gtf annotation

gtf_new <- gtf[gtf$transcript_id %in% trans_keep, ]


export(gtf_new, gtf_new_path)


######################## - kallisto expected counts


out_dir <- paste0(rwd, "2_counts/kallistofiltered5/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


counts_list <- lapply(1:length(samples), function(i){
  # i = 1
  print(i)
  
  abundance <- read.table(paste0(rwd, "2_counts/kallisto/", samples[i], "/", "abundance.txt"), header = TRUE, sep = "\t", as.is = TRUE)
  
  abundance <- abundance[abundance$target_id %in% trans_keep, ]
  
  counts <- data.frame(paste0(gene_trans[abundance[, "target_id"], "gene_id"], ":", abundance[, "target_id"]), counts = round(abundance[, "est_counts"]), stringsAsFactors = FALSE)

  colnames(counts) <- c("group_id", samples[i])
  
  write.table(counts, paste0(out_dir, samples[i], ".counts"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(NULL)
  
})


























