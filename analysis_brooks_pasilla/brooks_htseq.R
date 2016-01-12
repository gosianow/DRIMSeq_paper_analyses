######################################################
# BioC 2.14

# Created 13 Nov 2015 
# Get htseq counts

# Updated 24 Nov 2015
# Use kallisto reduced gtf

#######################################################

library("DEXSeq")
library("tools")

##############################################################################
# Test arguments
##############################################################################

# gtf='/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70_kallistoest_atleast5.gtf'
# count_method='htseqprefiltered5'

##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(gtf)
print(count_method)

##############################################################################

setwd("/home/Shared/data/seq/brooks_pasilla/")


DEXSeq_gff <- paste0(file_path_sans_ext(gtf), ".DEXSeq.flattened.rNO.gff")

counts_out <- paste0("2_counts/", count_method ,"/")
dir.create(counts_out, recursive = TRUE, showWarnings = FALSE)


##############################################################################
# metadata
##############################################################################


metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors=F, sep="\t", header=T) 

metadata


##############################################################################
# htseq counts
##############################################################################


# python scripts
pkgDir <- system.file(package="DEXSeq")
pythDir <- file.path(pkgDir, "python_scripts")
list.files(pythDir)


### crerate gff file # disable aggregation with the option “-r no”
# system(paste0("python ", pythDir, "/dexseq_prepare_annotation.py --help "))

python_cmd1 <- paste0("python ", pythDir, "/dexseq_prepare_annotation.py -r no ", gtf, " ", DEXSeq_gff)

cat(python_cmd1, fill = TRUE)

system(python_cmd1)



### exon counts # by name
# system(paste0("python ", pythDir, "/dexseq_count.py --help "))

python_cmd2 <- with(metadata, paste0("python ", pythDir, "/dexseq_count.py -p ", ifelse(LibraryLayout == "PAIRED", "yes", "no"), " -s no -f bam -r pos ", DEXSeq_gff, " 1_reads/tophat_2.0.14/", sampleName, "/accepted_hits.bam ", counts_out, sampleName , ".txt"))

cat(python_cmd2, fill = TRUE)




for(i in 1:nrow(metadata)){
  # i = 1 
  print(i)
  
  system(python_cmd2[i])
  
}



















