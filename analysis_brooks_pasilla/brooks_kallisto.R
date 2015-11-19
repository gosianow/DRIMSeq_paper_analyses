######################################################
# BioC 2.14

# Created 15 Nov 2015 

# Get kallisto counts

#######################################################

rwd <- "/home/Shared/data/seq/brooks_pasilla"

setwd(rwd)


######################################################################################################
# metadata
######################################################################################################

metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors=F, sep="\t", header=T) 

metadata


#########################################################################################
# Building an index
#########################################################################################

# kallisto index -i transcripts.idx transcripts.fasta.gz

cDNA.fasta <- "/home/Shared/data/annotation/Drosophila/Ensembl70/cDNA/Drosophila_melanogaster.BDGP5.70.cdna.all.fa"
index <- "/home/Shared/data/annotation/Drosophila/Ensembl70/kallisto/Drosophila_melanogaster.BDGP5.70.cdna.all.idx"

dir.create(dirname(index), recursive = TRUE, showWarnings = FALSE)


cmd <- paste("kallisto index -i",  index, cDNA.fasta, sep = " ")

system(cmd)



#########################################################################################
# Quantification
#########################################################################################
### Paired-end
# kallisto quant -i index -o output pairA_1.fastq pairA_2.fastq pairB_1.fastq pairB_2.fastq
### Single-end 
# kallisto quant -i index -o output --single -l 180 file1.fastq.gz file2.fastq.gz file3.fastq.gz

index <- "/home/Shared/data/annotation/Drosophila/Ensembl70/kallisto/Drosophila_melanogaster.BDGP5.70.cdna.all.idx"

fastq_dir <- paste0(rwd, "/1_reads/fastq/")


for(i in 3:nrow(metadata)){
  # i = 2
  
  sample <- metadata$SampleName[i]
  
  out_dir <- paste0(rwd, "/2_counts/kallisto/", sample, "/")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  switch(metadata$LibraryLayout[i], 
         
         SINGLE = {
           
           fastq <- paste0(fastq_dir, unlist(strsplit(metadata$fastq1[i], split = ",")), collapse = " ")
           
           cmd <- paste("kallisto quant -i", index, "-o", out_dir, "-b 0 -t 10 --single -l ", metadata$avgLength[i], " ", fastq)
           
           system(cmd)
           
         },
         
         PAIRED = {
           
           fastq <- paste0(fastq_dir, unlist(strsplit(metadata$fastq1[i], split = ",")), " " ,fastq_dir, unlist(strsplit(metadata$fastq2[i], split = ",")), collapse = " ")
           
           cmd <- paste("kallisto quant -i", index, "-o", out_dir, "-b 0 -t 10 ", fastq)
           
           system(cmd)
           
           
         }
  )
  
}





#########################################################################################
# Save abundance in files that can be an input for DEXSeq (estimated counts)
#########################################################################################

library(rtracklayer)

gtf_dir = "/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf"

gtf <- import(gtf_dir)

gt <- unique(mcols(gtf)[, c("gene_id", "transcript_id")])
rownames(gt) <- gt$transcript_id



############## save results in tables

out_dir <- paste0(rwd, "/2_counts/kallisto/")
samples <- metadata$SampleName

counts_list <- lapply(1:length(samples), function(i){
  # i = 1
  print(i)
  
  abundance <- read.table(paste0(out_dir, samples[i], "/", "abundance.txt"), header = TRUE, sep = "\t", as.is = TRUE)
  
  counts <- data.frame(paste0(gt[abundance$target_id, "gene_id"], ":", abundance$target_id), counts = round(abundance$est_counts), stringsAsFactors = FALSE)
  
  colnames(counts) <- c("group_id", samples[i])
  
  write.table(counts, paste0(out_dir, samples[i], ".counts"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(counts)
  
})



counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)


write.table(counts, paste0(out_dir, "kallisto_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)












































