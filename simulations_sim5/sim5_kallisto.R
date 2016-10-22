######################################################
## <<sim5_kallisto.R>>

# BioC 3.2
# Created 15 Nov 2015 
# Updated 18 Oct 2016

#######################################################
Sys.time()
#######################################################

library(tools)
library(rtracklayer)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/gosia/multinomial_project/simulations_sim5/hsapiens_withde_nonull'
# cDNA_fasta='0_annotation/TopHatTranscriptomeIndex/Homo_sapiens.GRCh37.71.dna.primary_assembly.fa'
# gtf_path='0_annotation/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding.gtf'
# conversion_path='0_annotation/KallistoIndex/TranscriptID_conversion.txt'
# out_dir='2_counts/kallisto'


##############################################################################
# Read in the arguments
##############################################################################

rm(list = ls())

args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

##############################################################################


setwd(rwd)

index <- paste0("0_annotation/KallistoIndex/", basename(file_path_sans_ext(cDNA_fasta)) ,".idx")

dir.create(dirname(index), recursive = TRUE, showWarnings = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

##############################################################################
# metadata
##############################################################################

metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata



##############################################################################
# Building an index
##############################################################################

# kallisto index -i transcripts.idx transcripts.fasta.gz

cmd <- paste("kallisto index -i",  index, cDNA_fasta, sep = " ")

system(cmd)



##############################################################################
# Quantification
##############################################################################
### Paired-end
# kallisto quant -i index -o output pairA_1.fastq pairA_2.fastq pairB_1.fastq pairB_2.fastq
### Single-end 
# kallisto quant -i index -o output --single -l 180 file1.fastq.gz file2.fastq.gz file3.fastq.gz


for(i in 1:nrow(metadata)){
  # i = 1

  out_dir_tmp <- paste0(rwd, "/", out_dir, "/sample", i, "/")
  dir.create(out_dir_tmp, recursive = TRUE, showWarnings = FALSE)

  fastq_dir <- paste0(rwd, "/1_reads/reads/")

  fastq <- paste0( fastq_dir, paste0("sample", i, "/sample_", i, "_1.fq "), fastq_dir, paste0("sample", i, "/sample_", i, "_2.fq"), collapse = " ")

  cmd <- paste("kallisto quant -i", index, "-o", out_dir_tmp, "-b 0 -t 5 ", fastq)

  system(cmd)

}



#########################################################################################
# Save the kallisto abundance (estimated counts) into files that can be used as an input for DEXSeq 
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
  
  abundance <- read.table(paste0(out_dir, "/sample", i, "/", "abundance.txt"), header = TRUE, sep = "\t", as.is = TRUE)
  
  abundance$target_id <- as.character(abundance$target_id)
  
  counts <- data.frame(paste0(id_conv[abundance$target_id, "gene_id"], ":", id_conv[abundance$target_id, "transcript_id"]), counts = round(abundance$est_counts), stringsAsFactors = FALSE)
  
  colnames(counts) <- c("group_id", paste0("sample_", i))
  
  write.table(counts, paste0(out_dir, "/kallisto", i, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(NULL)
  
})































##################################
# sim5_kallisto.R done!
##################################