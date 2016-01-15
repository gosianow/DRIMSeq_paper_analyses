######################################################
## ----- kim_kallisto
## <<kim_kallisto.R>>


# BioC 3.1
# Crated 19 June 2015
# Get kallisto counts

#######################################################

Sys.time()

#######################################################

library(tools)
library(rtracklayer)

##############################################################################
# Test arguments
##############################################################################


# rwd='/home/Shared/data/seq/kim_adenocarcinoma'
# cDNA_fasta='/home/Shared/data/annotation/Human/Ensembl_GRCh37.71/cDNA/Homo_sapiens.GRCh37.71.cdna.all.fa'
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

index <- paste0("0_annotation/kallisto/", basename(file_path_sans_ext(cDNA_fasta)) ,".idx")

dir.create(dirname(index), recursive = TRUE, showWarnings = FALSE)


##########################################################################
# load metadata
##########################################################################


metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors = FALSE, sep="\t", header=TRUE) 


sri <- read.table("3_metadata/SraRunInfo.csv", stringsAsFactors = F, sep = ",", header = T)

keep <- grep(" RNA-seq", sri$LibraryName)
sri <- sri[keep, ]



#########################################################################################
# Building an index
#########################################################################################

# kallisto index -i transcripts.idx transcripts.fasta.gz


cmd <- paste("kallisto index -i",  index, cDNA_fasta, sep = " ")

system(cmd)



#########################################################################################
# Quantification
#########################################################################################

# kallisto quant -i transcripts.idx -o output -b 100 -t 10 reads_1.fastq.gz reads_2.fastq.gz
# kallisto quant -i index -o output pairA_1.fastq pairA_2.fastq pairB_1.fastq pairB_2.fastq



fastq.path <- paste0(rwd, "1_reads/fastq/")

samples <- unique(sri$SampleName)


for(i in 1:length(samples)){
  # i = 1
  
  out.path <- paste0(rwd, "2_counts/kallisto/", samples[i], "/")
  dir.create(out.path, recursive = TRUE, showWarnings = FALSE)
  
  ssr <- sri[sri$SampleName == samples[i], "Run"]
  
  fastqList <- lapply(ssr, function(s){
    
    paste0(fastq.path, s,"_1.fastq.gz ", fastq.path, s,"_2.fastq.gz ")
    
  })
  
  fastq <- do.call(paste, fastqList)
  
  cmd <- paste("kallisto quant -i", index, "-o", out.path, "-b 0 -t 10", fastq) 
  
  system(cmd)
  
}





#########################################################################################
# Save abundance in files that can be an input for DEXSeq (estimated counts)
#########################################################################################


out.path.tmp <- paste0(rwd, "2_counts/kallisto/")

gtf <- import(gtf_path)

geneTrans <- unique(mcols(gtf)[, c("gene_id", "transcript_id")])
rownames(geneTrans) <- geneTrans$transcript_id



############## save results in tables

samples <- metadata$sampleName

counts_list <- lapply(1:length(samples), function(i){
  # i = 1
  print(i)
  
  abundance <- read.table(paste0(rwd, "2_counts/kallisto/", samples[i], "/", "abundance.txt"), header = TRUE, sep = "\t", as.is = TRUE)
  
  counts <- data.frame(paste0(geneTrans[abundance$target_id, "gene_id"], ":", abundance$target_id), counts = round(abundance$est_counts), stringsAsFactors = FALSE)
  
  colnames(counts) <- c("group_id", samples[i])
  
  write.table(counts, paste0(out.path.tmp, samples[i], ".counts"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(counts)
  
})
  


# counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)













































