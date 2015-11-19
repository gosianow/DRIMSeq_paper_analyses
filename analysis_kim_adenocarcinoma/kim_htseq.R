######################################################
# BioC 2.14
# Created 04 Nov 2014 

# Get htseq counts

#######################################################



setwd("/home/Shared/data/seq/kim_adenocarcinoma/")


library("DEXSeq")



######################################################################################################
# metadata
######################################################################################################

metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors = FALSE, sep="\t", header=TRUE) 



######################################################################################################
# htseq counts
######################################################################################################

counts_out <- "2_counts/htseq/"
dir.create(counts_out, recursive = TRUE)


gtf= "/home/Shared_penticton/data/annotation/Human/Ensembl_GRCh37.71/gtf/Homo_sapiens.GRCh37.71.gtf"

DEXSeq.gff = "/home/Shared/data/annotation/Human/Ensembl_GRCh37.71/DEXSeq_1.10.8_gff/Homo_sapiens.GRCh37.71.DEXSeq.flattened.rNO.gff"

# python scripts
pkgDir <- system.file(package="DEXSeq")
pythDir <- file.path(pkgDir, "python_scripts")
list.files(pythDir)


### crerate gff file # disable aggregation with the option “-r no”
system(paste0("python ", pythDir, "/dexseq_prepare_annotation.py --help "))

python.cmd1 <- paste0("python ", pythDir, "/dexseq_prepare_annotation.py -r no ", gtf, " ", DEXSeq.gff)
cat(python.cmd1)
system(python.cmd1)



### add chr to gff
## awk '{print "chr"$0}' Homo_sapiens.GRCh37.71.DEXSeq.flattened.rNO.gff > Homo_sapiens.GRCh37.71.DEXSeq.flattened.rNO.chr.gff

DEXSeq.gff <- "/home/Shared/data/annotation/Human/Ensembl_GRCh37.71/DEXSeq_1.10.8_gff/Homo_sapiens.GRCh37.71.DEXSeq.flattened.rNO.chr.gff"


### exon counts # by name
system(paste0("python ", pythDir, "/dexseq_count.py --help "))

python.cmd2 <- with(metadata, paste0("python ", pythDir, "/dexseq_count.py -p yes -s no -f bam -r pos ", DEXSeq.gff, " bam_insilicodb/", ids, "_s.bam ", counts_out, ids , ".counts \n"))
cat(python.cmd2)

for(i in 1:nrow(metadata)){
  # i = 1 
  system(python.cmd2[i])
  
}



### Merge counts into one file

counts_list <- lapply(1:nrow(metadata), function(i){
  # i = 1
  
  cts <- read.table(paste0(counts_out, metadata$ids[i], ".counts"), header = FALSE, as.is = TRUE)
  colnames(cts) <- c("group_id", metadata$ids[i])
  
  return(cts)
  
})


counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)

tail(counts)
counts <- counts[!grepl(pattern = "_", counts$group_id), ]


write.table(counts, paste0(counts_out, "htseq_counts.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)











