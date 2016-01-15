######################################################
## ----- brooks_download
## <<brooks_download.R>>


# BioC 2.14
# Created 13 Nov 2015

#######################################################

Sys.time()

#######################################################



##############################################################################
# Test arguments
##############################################################################


# rwd='/home/Shared/data/seq/brooks_pasilla'
# gtf='/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf'
# DNA_fasta='/home/Shared/data/annotation/Drosophila/Ensembl70/genome/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa'


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

stopifnot(file.exists("3_metadata/SraRunInfo.csv"))

bowtie_ind_dir <- "0_annotation/bowtie2_index"
dir.create(bowtie_ind_dir, recursive = TRUE)


##############################################################################
# create metadata file -  metadata
##############################################################################


fastq_dir <- paste0("1_reads/fastq/")
dir.create(fastq_dir, recursive = TRUE)


sri.org <- read.table("3_metadata/SraRunInfo.csv", stringsAsFactors=F, sep=",", header=T)
keep <- grep(paste("GSM4611", 76:82, sep="", collapse="|"), sri.org$SampleName)
sri <- sri.org[keep,]

sri$LibraryName = gsub("S2_DRSC_","",sri$LibraryName)
metadata = unique(sri[,c("LibraryName","LibraryLayout", "SampleName")])

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


write.table(metadata, "3_metadata/metadata.xls", sep="\t", row.names=F, quote=F)





##############################################################################
# download data
##############################################################################


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
# Bowtie2 index
####################################


cmd <- paste0("bowtie2-build -f ", DNA_fasta, " ", bowtie_ind_dir, "/Dme1_BDGP5_70")
system(cmd)


bowtie_ind <- paste0(bowtie_ind_dir, "/Dme1_BDGP5_70")


####################################
# TopHat2
####################################


dir.create("1_reads/tophat_2.0.14", recursive = TRUE)

cmd <- with(metadata, paste("tophat2 -G", gtf, "-p 10 -o", paste0("1_reads/tophat_2.0.14/", SampleName), bowtie_ind, fastq1, fastq2, "\n"))

cat(cmd)

for(i in 1:length(cmd))
  system.time(system(cmd[i]))



####################################
# organize BAM files with samtools
####################################

# dir.create("1_reads/sam/", recursive = TRUE)

# for (i in 1:nrow(metadata)){
#   # i = 2
#   sam.in <- paste0("1_reads/tophat_2.0.14/" ,metadata$SampleName[i], "/accepted_hits.bam")
#   sam.out <- paste0("1_reads/sam/" ,metadata$SampleName[i])
  
#   #   sort by position and index for IGV and htseq
#   cmdt = paste0("samtools sort ", sam.in, " ", sam.out, "_s", "\n",
#                 "samtools index ", sam.out, "_s.bam", "\n",
#                 "samtools view -o ", sam.out, "_s.sam ", sam.out, "_s.bam", "\n")
  
#   cat(cmdt)
#   system(cmdt)
# }



