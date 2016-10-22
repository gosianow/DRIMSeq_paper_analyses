######################################################
## <<sim5_kallisto_filter_gtf_fa.R>>

# BioC 3.2
# Crated 24 Nov 2015
# Updated 18 Oct 2016

# Use TPMs fractions to filter out the low expressed transcripts from, gtf, fa and kallisto counts 

#######################################################
Sys.time()
#######################################################

library(DRIMSeq) # used to do the filtering
library(limma)
library(rtracklayer)
library(tools)
library(Biostrings) # for fasta

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/gosia/multinomial_project/simulations_sim5/hsapiens_withde_nonull'
# fa_path='0_annotation/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa'
# gtf_path='0_annotation/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf'
# conversion_path='0_annotation/KallistoIndex/TranscriptID_conversion.txt'
# kallisto_dir='2_counts/kallisto'
# out_dir='2_counts/kallisto_txfilt_5'

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

gtf_new_path <- paste0("0_annotation/INCOMPLETE_KALLISTOEST/", basename(file_path_sans_ext(gtf_path)) ,"_kallistoest_atleast5.gtf")
fa_new_path <- paste0("0_annotation/INCOMPLETE_KALLISTOEST/", basename(file_path_sans_ext(fa_path)) ,"_kallistoest_atleast5.fa")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

##########################################################################
# load metadata
##########################################################################

metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata



#########################################################################################
### Load kallisto TPMs
#########################################################################################

gtf <- import(gtf_path)

id_conv <- read.table(conversion_path, header = FALSE, sep = " ", as.is = TRUE)
colnames(id_conv) <- c("target_id", "transcript_id")
rownames(id_conv) <- id_conv[, "target_id"]

# add gene ids
gt <- unique(mcols(gtf)[, c("gene_id", "transcript_id")])

all(id_conv$transcript_id %in% gt$transcript_id)

mm <- match(id_conv$transcript_id, gt$transcript_id)
id_conv$gene_id <- NA
id_conv$gene_id[!is.na(mm)] <- gt[mm[!is.na(mm)], "gene_id"]


### load TPMs
counts_list <- lapply(1:nrow(metadata), function(i){
  # i = 1
  
  abundance <- read.table(paste0(kallisto_dir, "/sample", i, "/", "abundance.txt"), header = TRUE, sep = "\t", as.is = TRUE)
  abundance$target_id <- as.character(abundance$target_id)
  
  counts <- data.frame(paste0(id_conv[abundance$target_id, "gene_id"], ":", id_conv[abundance$target_id, "transcript_id"]), counts = abundance$tpm, stringsAsFactors = FALSE)
  
  colnames(counts) <- c("group_id", paste0("sample_", i))
  
  return(counts)
  
})


tpm <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)

group_split <- strsplit2(tpm[, "group_id"], ":")

d_org <- dmDSdata(counts = tpm[, -1], gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = colnames(tpm[, -1]), group = rep("C1", ncol(tpm[, -1])))


d_filt <- dmFilter(d_org, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 1, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0.05)

trans_keep <- counts(d_filt)$feature_id


#########################################################################################
### Subset gtf
#########################################################################################

gtf_new <- gtf[gtf$transcript_id %in% trans_keep, ]

export(gtf_new, gtf_new_path)


###############################################################################
### Subset fasta
###############################################################################


fa <- readDNAStringSet(fa_path)

fa_names <- sapply(strsplit(names(fa), " "), .subset, 2)

mm <- match(trans_keep, fa_names)

fa_new <- fa[mm]

writeXStringSet(fa_new, filepath = fa_new_path)



#########################################################################################
### Subset kallisto expected counts
#########################################################################################


counts_list <- lapply(1:nrow(metadata), function(i){
  # i = 1

  abundance <- read.table(paste0(kallisto_dir, "/kallisto", i, ".txt"), header = FALSE, sep = "\t", as.is = TRUE)
  
  trans_name <- strsplit2(abundance[, 1], ":")[, 2]
  
  abundance <- abundance[trans_name %in% trans_keep, ]
  
  write.table(abundance, paste0(out_dir, "/kallisto", i, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(NULL)
  
})















