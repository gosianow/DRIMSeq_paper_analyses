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


fastq_dir <- paste0("1_reads/fastq/")

sri.org <- read.table("3_metadata/SraRunInfo.csv", stringsAsFactors=F, sep=",", header=T)
keep <- grep(paste("GSM4611", 76:82, sep="", collapse="|"), sri.org$SampleName)
sri <- sri.org[keep,]

sri$LibraryName = gsub("S2_DRSC_","",sri$LibraryName)
metadata = unique(sri[,c("LibraryName","LibraryLayout", "SampleName", "avgLength")])

for(i in seq_len(nrow(metadata))) {
  rw = (sri$LibraryName==metadata$LibraryName[i])
  if(metadata$LibraryLayout[i]=="PAIRED") {
    metadata$fastq1[i] = paste0(fastq_dir, sri$Run[rw],"_1.fastq.gz",collapse=",")
    metadata$fastq2[i] = paste0(fastq_dir, sri$Run[rw],"_2.fastq.gz",collapse=",")
    # metadata$ReadLength[i] <- metadata$avgLength[i] / 2
    # metadata$SRR[i] <- paste0(sri$Run[rw], collapse=",")
  } else {
    metadata$fastq1[i] = paste0(fastq_dir, sri$Run[rw],".fastq.gz",collapse=",")
    metadata$fastq2[i] = ""
    # metadata$ReadLength[i] <- metadata$avgLength[i]
    # metadata$SRR[i] <- paste0(sri$Run[rw], collapse=",")
  }
}

metadata$condition = "CTL"
metadata$condition[grep("CG8144_RNAi",metadata$LibraryName)] = "KD"
metadata$shortname = paste( seq_len(nrow(metadata)), substr(metadata$condition,1,2),  substr(metadata$LibraryLayout,1,2), metadata$ReadLength, sep=".")
metadata$sampleName <- metadata$SampleName
metadata$UniqueName <- paste0(1:nrow(metadata), "_", metadata$SampleName)

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

fastq_dir <- paste0(rwd, "/")


for(i in 1:nrow(metadata)){
  # i = 2
  
  sample <- metadata$UniqueName[i]
  
  out_dir <- paste0(rwd, "/2_counts/kallisto/", sample, "/")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  switch(metadata$LibraryLayout[i], 
         
         SINGLE = {
           
           fastq <- paste0(fastq_dir, unlist(strsplit(metadata$fastq1[i], split = ",")), collapse = " ")
           
           cmd <- paste("kallisto quant -i", index, "-o", out_dir, "-b 0 -t 5 --single -l ", metadata$avgLength[i], fastq)
           
           system(cmd)
           
         },
         
         PAIRED = {
           
           fastq <- paste0(fastq_dir, unlist(strsplit(metadata$fastq1[i], split = ",")), " " ,fastq_dir, unlist(strsplit(metadata$fastq2[i], split = ",")), collapse = " ")
           
           cmd <- paste("kallisto quant -i", index, "-o", out_dir, "-b 0 -t 5 ", fastq)
           
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
samples <- unique(metadata$SampleName)


counts_list <- lapply(1:length(samples), function(i){
  # i = 4
  print(i)
  
  indx <- which(metadata$SampleName == samples[i])
  
  if(length(indx) == 1){
    
    abundance <- read.table(paste0(out_dir, metadata$UniqueName[indx], "/", "abundance.txt"), header = TRUE, sep = "\t", as.is = TRUE)
    
  }else{
    
    abundance <- lapply(indx, function(j){
      # j = indx[1]
      
      abundance_tmp <- read.table(paste0(out_dir, metadata$UniqueName[j], "/", "abundance.txt"), header = TRUE, sep = "\t", as.is = TRUE)
      abundance_tmp <- abundance_tmp[, c("target_id", "est_counts")]
      abundance_tmp
      
      })
    
    abundance <- Reduce(function(...) merge(..., by = "target_id", all = TRUE, sort = FALSE), abundance)
    est_counts <- rowSums(abundance[, -1])
    
    abundance <- data.frame(target_id = abundance$target_id, est_counts = est_counts, stringsAsFactors = FALSE)
    
  }

  counts <- data.frame(paste0(gt[abundance$target_id, "gene_id"], ":", abundance$target_id), counts = round(abundance$est_counts), stringsAsFactors = FALSE)
  
  colnames(counts) <- c("group_id", samples[i])
  
  write.table(counts, paste0(out_dir, samples[i], ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(counts)
  
})



# counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)
# write.table(counts, paste0(out_dir, "kallisto_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)












































