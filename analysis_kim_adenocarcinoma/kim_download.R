######################################################
## ----- kim_download
## <<kim_download.R>>

# BioC 2.14
# Created 15 Jan 2015 

##############################################################################

Sys.time()

##############################################################################


##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/kim_adenocarcinoma'

##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

print(rwd)


##############################################################################

setwd(rwd)


##############################################################################
# prepare metadata file
##############################################################################


metadata <- read.table("3_metadata/Malgorzata Nowicka2014-11-04GSE37764.csv", stringsAsFactors=F, sep=",", header=T) 
metadata <- metadata[metadata$X == "RNA-seq",]

metadata$sampleName <- metadata$ids
metadata$condition <- metadata$Tissue.Type

write.table(metadata, "3_metadata/metadata.xls", sep="\t", row.names=F, quote=F)


##############################################################################
# download bam files from insilicodb
##############################################################################



#### download tophat results from insilicodb ####


link <- paste0("https://insilicodb.com/app/publicutilities/getfile?seriesName=GSE37764&platformName=GPL10999&file=source/GSE37764-measurements/", metadata$ids ,"/rnaseq/v0.9/tophat_out/accepted_hits.bam")

out_dir <- paste0("1_reads/tophat_insilicodb/", metadata$ids)



for(i in 1:nrow(metadata)){
  print(i)
  
  dir.create(out_dir[i], recursive = TRUE)
  
  download.file(link[i], paste0(out_dir[i], "/accepted_hits.bam"), method="wget")
  
}



#### organize BAM files with samtools  - not needed anymore

# dir.create(paste0("1_reads/bam_insilicodb/"), recursive = TRUE)
# 
# for (i in 1:nrow(metadata)){
#   # i=1
#   bam.in <- paste0("1_reads/tophat_insilicodb/" , metadata$ids[i], "/accepted_hits.bam")
#   bam.out <- paste0("1_reads/bam_insilicodb/" ,metadata$ids[i])
#   
#   #   sort by position (for featureCounts)
#   #   convert to SAM (for DEXSeq) not needed 
#   cmdt = paste0("samtools sort ", bam.in, " ", bam.out, "_s", "\n",
#                 "samtools view -o ", bam.out, "_s.sam ", bam.out, "_s.bam", "\n" )
#   
#   #   sort by name (for DEXSeq) not needed
#   cmdt <- paste0("samtools sort -n ", bam.in, " ", bam.out, "_sn", "\n")
# 
#   #   add index for IGV
#   cmdt <- paste0("samtools index ", bam.out, "_s.bam", "\n")
#   
#   
#   cat(cmdt)
#   system(cmdt)
#   
# }
# 



##############################################################################
# download sra files
##############################################################################


out_dir <- "1_reads/fastq/"
dir.create(out_dir)


sri <- read.table("3_metadata/SraRunInfo.csv", stringsAsFactors = F, sep = ",", header = T)

keep <- grep(" RNA-seq", sri$LibraryName)
sri <- sri[keep, ]

# "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR493/SRR493937/SRR493937.sra"

download_path <- paste0("ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/", substr(sri$Run, 1, 6), "/", sri$Run, "/", sri$Run, ".sra" )


files.sra <- paste(out_dir, basename(download_path), sep = "")


for (i in 1:length(files.sra)) 
  download.file(download_path[i], files.sra[i])



#### convert sra to fastq

for (i in 1:length(files.sra)) {
  cmd <- paste0("fastq-dump -O ", out_dir, " --split-3 ", files.sra[i])
  cat(cmd, "\n")
  system(cmd)
}



system(paste0("gzip ", out_dir, "*.fastq"))
system(paste0("/bin/rm ", out_dir, "*.sra"))
















































