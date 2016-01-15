######################################################
## ----- brooks_kallisto_filter_gtf
## <<brooks_kallisto_filter_gtf.R>>


# BioC 3.2

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


# rwd='/home/Shared/data/seq/brooks_pasilla'
# gtf_path='/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf'



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
# Use TPMs fractions to filter out the low expressed transcripts from 
# - gtf annotation
# - kallisto expected counts
#########################################################################################


gtf <- import(gtf_path)

gt <- unique(mcols(gtf)[, c("gene_id", "transcript_id")])
rownames(gt) <- gt$transcript_id


### load TPMs


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
      abundance_tmp <- abundance_tmp[, c("target_id", "tpm")]
      abundance_tmp
      
    })
    
    abundance <- Reduce(function(...) merge(..., by = "target_id", all = TRUE, sort = FALSE), abundance)
    tpm <- rowSums(abundance[, -1])
    
    abundance <- data.frame(target_id = abundance$target_id, tpm = tpm, stringsAsFactors = FALSE)
    
  }
  
  counts <- data.frame(paste0(gt[abundance$target_id, "gene_id"], ":", abundance$target_id), counts = abundance$tpm, stringsAsFactors = FALSE)
  
  colnames(counts) <- c("group_id", samples[i])
  
  return(counts)
  
})


tpm <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)

group_split <- strsplit2(tpm[, "group_id"], ":")

d_org <- dmDSdata(counts = tpm[, -1], gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = colnames(tpm[, -1]), group = rep("C1", ncol(tpm[, -1])))


d_filt <- dmFilter(d_org, min_samps_gene_expr = 0, min_samps_feature_expr = 0, min_samps_feature_prop = 1, min_gene_expr = 0, min_feature_expr = 0, min_feature_prop = 0.05)

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
  
  abundance <- read.table(paste0(rwd, "2_counts/kallisto/", samples[i], ".txt"), header = FALSE, sep = "\t", as.is = TRUE)
  
  trans <- strsplit2(abundance[, 1], ":")[, 2]
  
  abundance <- abundance[trans %in% trans_keep, ]
  
  write.table(abundance, paste0(out_dir, samples[i], ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(NULL)
  
})


























