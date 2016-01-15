##############################################################################
# BioC 2.14

# Created 04 Nov 2014 
# Get htseq counts

##############################################################################

library(DEXSeq)
library(tools)
library(BiocParallel)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/kim_adenocarcinoma'
# gtf='/home/Shared_penticton/data/annotation/Human/Ensembl_GRCh37.71/gtf/Homo_sapiens.GRCh37.71_kallistoest_atleast5.gtf'
# count_method='htseqprefiltered5'

##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

print(gtf)
print(count_method)

##############################################################################

setwd(rwd)


DEXSeq_gff <- paste0(file_path_sans_ext(gtf), ".DEXSeq.flattened.rNO.gff")
DEXSeq_gff_chr <- paste0(file_path_sans_ext(gtf), ".DEXSeq.flattened.rNO.chr.gff")

counts_out <- paste0("2_counts/", count_method ,"/")
dir.create(counts_out, recursive = TRUE, showWarnings = FALSE)


#############################################################################
# metadata
#############################################################################

metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors = FALSE, sep="\t", header=TRUE) 
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
cat(python_cmd1)

system(python_cmd1)



### add chr to gff
cmd <- paste0("awk '{print \"chr\"$0}' ", DEXSeq_gff ," > ", DEXSeq_gff_chr)

cat(cmd, fill = TRUE)

system(cmd)



### exon counts # by name
# system(paste0("python ", pythDir, "/dexseq_count.py --help "))

python_cmd2 <- with(metadata, paste0("python ", pythDir, "/dexseq_count.py -p yes -s no -f bam -r pos ", DEXSeq_gff_chr, " 1_reads/tophat_insilicodb/", sampleName, "/accepted_hits.bam ", counts_out, sampleName , ".counts \n"))

cat(python_cmd2, fill = TRUE)



# for(i in 1:nrow(metadata)){
#   # i = 1 
#   system(python_cmd2[i])
#   
# }


bplapply(1:nrow(metadata), function(i){
  
  system(python_cmd2[i])
  
  return(NULL)
  
}, BPPARAM = MulticoreParam(workers = 4))


# ### Merge counts into one file
# 
# counts_list <- lapply(1:nrow(metadata), function(i){
#   # i = 1
#   
#   cts <- read.table(paste0(counts_out, metadata$sampleName[i], ".counts"), header = FALSE, as.is = TRUE)
#   colnames(cts) <- c("group_id", metadata$sampleName[i])
#   
#   return(cts)
#   
# })
# 
# 
# counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)
# 
# tail(counts)
# counts <- counts[!grepl(pattern = "_", counts$group_id), ]
# 
# 
# write.table(counts, paste0(counts_out, "htseq_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)











