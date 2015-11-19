# BioC 3.1
# Created 13 Nov 2015

##############################################################################
# create metadata file -  metadata
##############################################################################

setwd("/home/Shared/data/seq/brooks_pasilla/")



sri.org <- read.table("3_metadata/SraRunInfo.csv", stringsAsFactors=F, sep=",", header=T)
keep <- grep(paste("GSM4611", 76:82, sep="", collapse="|"), sri.org$SampleName)
sri <- sri.org[keep,]

# manual trimming of reads to have equal lenght 
which.trim <- sri[sri$SampleName=="GSM461179", c("Run")]

# cmd <- paste0("java -jar /usr/local/software/Trimmomatic-0.30/trimmomatic-0.30.jar  SE -threads 20 download/", which.trim , ".fastq.gz download/", which.trim , ".trimed.fastq.gz CROP:40 \n")
# cat(cmd)
# for(i in 1:length(cmd))
#   system(cmd[i])


sri[sri$Run %in% which.trim, "avgLength"] <- 40
sri[sri$Run %in% which.trim, "Run"] <- paste0(sri[sri$Run %in% which.trim, "Run"], ".trimed")

sri$LibraryName = gsub("S2_DRSC_","",sri$LibraryName) # trim label
metadata = unique(sri[,c("LibraryName","LibraryLayout", "SampleName", "avgLength" )])

for(i in seq_len(nrow(metadata))) {
  rw = (sri$LibraryName==metadata$LibraryName[i])
  if(metadata$LibraryLayout[i]=="PAIRED") {
    metadata$fastq1[i] = paste0(sri$Run[rw],"_1.fastq.gz",collapse=",")
    metadata$fastq2[i] = paste0(sri$Run[rw],"_2.fastq.gz",collapse=",")
    metadata$ReadLength[i] <- metadata$avgLength[i] / 2
    metadata$SRR[i] <- paste0(sri$Run[rw], collapse=",")
  } else {
    metadata$fastq1[i] = paste0(sri$Run[rw],".fastq.gz",collapse=",")
    metadata$fastq2[i] = ""
    metadata$ReadLength[i] <- metadata$avgLength[i]
    metadata$SRR[i] <- paste0(sri$Run[rw], collapse=",")
  }
}

metadata$condition = "CTL"
metadata$condition[grep("RNAi",metadata$LibraryName)] = "KD"
metadata$shortname = paste( seq_len(nrow(metadata)), substr(metadata$condition,1,2),  substr(metadata$LibraryLayout,1,2), metadata$ReadLength, sep=".")
metadata$sampleName <- metadata$SampleName


write.table(metadata, "3_metadata/metadata.xls", sep="\t", row.names=F, quote=F)





##############################################################################
# download data
##############################################################################

setwd("/home/Shared/data/seq/brooks_pasilla/")

sri <- read.table("3_metadata/SraRunInfo.csv", stringsAsFactors=F, sep=",", header=T)
keep <- grep(paste("GSM4611", 76:82, sep="", collapse="|"), sri$SampleName)
sri <- sri[keep,]


files.sra <- paste("1_reads/fastq/", basename(sri$download_path), sep="") 

for(i in 1:nrow(sri))
  download.file(sri$download_path[i], files.sra[i])


for(i in 1:nrow(sri)) {
  cmd <- paste0("fastq-dump -O 1_reads/fastq --split-3 ", files.sra[i])
  cat(cmd,"\n")
  system(cmd)
}

system("gzip 1_reads/fastq/*.fastq")
system("/bin/rm 1_reads/fastq/*.sra")



####################################
# TopHat2
####################################

gtf="/home/Shared_penticton/data/annotation/Drosophila/gtf_Ensembl70/Drosophila_melanogaster.BDGP5.70.gtf"
bow.ind="/home/Shared_penticton/data/annotation/Drosophila/bowtie2_index_Ensembl70/Dme1_BDGP5_70"

dir.create("1_reads/tophat_2.0.9", recursive = TRUE)

cmd = with(metadata, paste("tophat2 -G", gtf, "-p 20 -o", paste0("1_reads/tophat_2.0.9/", SampleName), bow.ind, fastq1, fastq2, "\n"))

cat(cmd)

for(i in 1:length(cmd))
  system.time(system(cmd[i]))



####################################
# organize BAM files with samtools
####################################


for (i in 1:nrow(metadata)){
  # i=2
  sam.in <- paste0("1_reads/tophat_2.0.9/" ,metadata$SampleName[i], "/accepted_hits.bam")
  sam.out <- paste0("1_reads/sam/" ,metadata$SampleName[i])
  
  #   sort by position and index for IGV and htseq
  cmdt = paste0("samtools sort ", sam.in, " ", sam.out, "_s", "\n",
                "samtools index ", sam.out, "_s.bam", "\n",
                "samtools view -o ", sam.out, "_s.sam ", sam.out, "_s.bam", "\n")
  
  cat(cmdt)
  system(cmdt)
}



