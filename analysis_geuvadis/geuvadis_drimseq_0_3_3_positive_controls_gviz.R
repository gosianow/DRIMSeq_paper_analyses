######################################################
## <<geuvadis_drimseq_0_3_3_positive_controls_gviz.R>>

# BioC 3.2
# Created 2 Mar 2016
# Modified 27 Apr 2016

# Plot gene structure for validated sQTLs

##############################################################################

Sys.time()

##############################################################################

library(DRIMSeq)
library(ggplot2)
library(reshape2)
library(iCOBRA)
library(plyr)
library(Gviz)
library(GenomicFeatures)
library(tools)
library(rtracklayer)
library(GenomicRanges)
library(limma)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/geuvadis'
# population='CEU'
# path_gtf='geuvadis_annotation/gencode.v12.annotation.gtf'
# valid_path='data/validation/glimmps/glimmps_valid_pcr.txt'
# positive_controls_out='drimseq_0_3_3_positive_controls_permutations_all_genes'

##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}

print(args)

##############################################################################

setwd(rwd)


positive_controls_out <- paste0(positive_controls_out, "/", basename(file_path_sans_ext(valid_path)), "/")

dir.create(positive_controls_out, recursive = TRUE, showWarnings = FALSE)

out_dir <- positive_controls_out

dir.create(paste0(out_dir, "figures/"), recursive = TRUE, showWarnings = FALSE)


##############################################################################
# validated genes
##############################################################################

valid <- read.table(valid_path, header = TRUE, sep = "\t", as.is = TRUE) 

valid 


##########################################################################
### Gene structure plot for validated snps
##########################################################################



gat <- GenomeAxisTrack()

gtf <- import(path_gtf)

txdb <- makeTxDbFromGFF(path_gtf, format="gtf")
genes <- genes(txdb)


txTr_global <- GeneRegionTrack(txdb, name = "transcripts",  transcriptAnnotation = "transcript", just.group = "above", cex.group = 1)


for(j in 1:nrow(valid)){
  # j = 2
  
  g <- valid$gene_id[j]
  chromosome <- as.character(seqnames(genes[g, ]))
  astart <- start(genes[g, ]) - 5000
  aend <- end(genes[g, ]) + 5000
  
  # ### Use gtf to define the gene structure
  # gtf_sub <- gtf[mcols(gtf)$gene_id == g & mcols(gtf)$type != "transcript" & mcols(gtf)$type != "gene", ]
  # mcols(gtf_sub)$gene <- mcols(gtf_sub)$gene_id
  # mcols(gtf_sub)$transcript <- mcols(gtf_sub)$transcript_id
  # mcols(gtf_sub)$feature <- as.character(mcols(gtf_sub)$type)
  # 
  # txTr <- GeneRegionTrack(gtf_sub, name = "transcripts",  transcriptAnnotation = "transcript", just.group = "above", cex.group = 1)
  # 
  # 
  # ht <- HighlightTrack(trackList = list(gat, txTr), start = c(valid$snp_position[j]), width = 1, chromosome = chromosome)
  # 
  # 
  # pdf(paste0(out_dir, "figures/gviz_", j, "_", valid$gene_name[j], ".pdf"), 12, 8)
  # 
  # plotTracks(c(ht), from = astart, to = aend, chromosome = chromosome)
  # 
  # dev.off()
  # 
  # 
  # 
  # pdf(paste0(out_dir, "figures/gviz_", j, "_", valid$gene_name[j], "_zoom.pdf"), 12, 8)
  # 
  # plotTracks(c(ht), from = valid$snp_position[j] - 5000, to = valid$snp_position[j] + 5000, chromosome = chromosome)
  # 
  # dev.off()
  
  if(!(valid$snp_position[j] < aend && valid$snp_position[j] > astart))
    next
  
  
  if(all(c("target_exon_start", "target_exon_end") %in% colnames(valid))){
    
    ht <- HighlightTrack(trackList = list(gat, txTr_global), start = c(valid$target_exon_start[j], valid$snp_position[j]), end = c(valid$target_exon_end[j], valid$snp_position[j] + 1), chromosome = chromosome, col = c("red", "blue"), lty = c(1, 2), lwd = c(2, 2))
    
  }else{

    ht <- HighlightTrack(trackList = list(gat, txTr_global), start = c(valid$snp_position[j]), end = c(valid$snp_position[j] + 1), chromosome = chromosome, col = "blue", lty = 1, lwd = 2)
    
  }
  
  
  pdf(paste0(out_dir, "figures/gviz_", j, "_", valid$gene_name[j], ".pdf"), 12, 8)
  
  plotTracks(c(ht), from = astart, to = aend, chromosome = chromosome)
  
  dev.off()
  
  
  
  pdf(paste0(out_dir, "figures/gviz_", j, "_", valid$gene_name[j], "_zoom.pdf"), 12, 8)
  
  plotTracks(c(ht), from = valid$snp_position[j] - 3000, to = valid$snp_position[j] + 3000, chromosome = chromosome)
  
  dev.off()
  
  
}







# data(cyp2b10)
# 
# 
# ## Construct the object
# grTrack <- GeneRegionTrack(start=26682683, end=26711643,
#   rstart=cyp2b10$start, rends=cyp2b10$end, chromosome=7, genome="mm9",
#   transcript=cyp2b10$transcript, gene=cyp2b10$gene, symbol=cyp2b10$symbol,
#   feature=cyp2b10$feature, exon=cyp2b10$exon,
#   name="Cyp2b10", strand=cyp2b10$strand)
# 
# ## Directly from the data.frame
# grTrack <- GeneRegionTrack(cyp2b10)












sessionInfo()

