######################################################
## <<sim5_htseq.R>>

# BioC 2.14

# Created 13 Nov 2015 
# Updated 18 Oct 2016

#######################################################

library("DEXSeq")
library("tools")

##############################################################################
# Test arguments
##############################################################################


# rwd='/home/gosia/multinomial_project/simulations_sim5/hsapiens_withde_nonull'
# gtf_path='0_annotation/INCOMPLETE_KALLISTOEST/Homo_sapiens.GRCh37.71.primary_assembly.protein_coding_kallistoest_atleast5.gtf'
# out_dir='2_counts/INCOMPLETE_KALLISTOEST/dexseq_nomerge_kallistoest_atleast5'


##############################################################################
# Read in the arguments
##############################################################################

rm(list = ls())

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

##############################################################################

setwd(rwd)


DEXSeq_gff <- paste0(file_path_sans_ext(gtf_path), ".flattened.nomerge.gff")
DEXSeq_gff

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


##############################################################################
# metadata
##############################################################################

metadata <- data.frame(sample_id = paste0("sample_",1:6), group=c(rep("C1", 3), rep("C2", 3)), stringsAsFactors = FALSE)
metadata$group <- as.factor(metadata$group)

metadata


##############################################################################
# htseq counts
##############################################################################


# python scripts
pkgDir <- system.file(package="DEXSeq")
pythDir <- file.path(pkgDir, "python_scripts")
list.files(pythDir)


### crerate gff file and disable aggregation with the option “-r no”

# system(paste0("python ", pythDir, "/dexseq_prepare_annotation.py --help "))

python_cmd1 <- paste0("python ", pythDir, "/dexseq_prepare_annotation.py -r no ", gtf_path, " ", DEXSeq_gff)

cat(python_cmd1, fill = TRUE)

system(python_cmd1)



### exon counts
# system(paste0("python ", pythDir, "/dexseq_count.py --help "))

python_cmd2 <- paste0("python ", pythDir, "/dexseq_count.py -p yes -s no -f bam -r pos ", DEXSeq_gff, " 1_reads/tophat/sample", 1:6, "/accepted_hits.bam ", out_dir, "/dexseq", 1:6, ".txt")

cat(python_cmd2, fill = TRUE)



for(i in 1:nrow(metadata)){
  # i = 1 
  print(i)
  
  system(python_cmd2[i])
  
}



















