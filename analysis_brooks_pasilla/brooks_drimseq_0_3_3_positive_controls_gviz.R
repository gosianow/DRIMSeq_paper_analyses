######################################################
## <<brooks_drimseq_0_3_3_positive_controls_gviz.R>>

# BioC 3.2
# Created 15 Jan 2015 
# Modified 27 Apr 2016

##############################################################################

Sys.time()

##############################################################################

library(plyr)
library(ggplot2)
library(rtracklayer)
library(reshape2)
library(Gviz)
library(GenomicFeatures)
library(tools)
library(limma)

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/brooks_pasilla'
# path_gtf='/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf'
# path_gtf_filtered = '/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70_kallistoest_atleast5.gtf'
# valid_path='5_validation/brooks_validated_genes.txt'
# method_out='drimseq_0_3_3'
# comparison_out='drimseq_0_3_3_positive_controls'

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

comparison_out <- paste0(comparison_out, "/")

dir.create(comparison_out, showWarnings = FALSE, recursive = TRUE)

dir.create(paste0(comparison_out, "figures/"), showWarnings = FALSE, recursive = TRUE)

path_gtf_dexseq <- paste0(file_path_sans_ext(path_gtf), ".DEXSeq.flattened.rNO.gff")
path_gtf_filtered_dexseq <- paste0(file_path_sans_ext(path_gtf_filtered), ".DEXSeq.flattened.rNO.gff")


##############################################################################
# metadata
##############################################################################

metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors=F, sep="\t", header=T) 

metadata

##############################################################################
# validated genes
##############################################################################

valid <- read.table(valid_path, header = TRUE, sep = "\t", as.is = TRUE) 

valid 


##############################################################################
# Plot coverage and annotations for the validated genes
##############################################################################

dir.create(paste0(comparison_out, "gviz/"))


options(ucscChromosomeNames=FALSE)

bam_dir <- "1_reads/tophat_2.0.14/"

metadata <- metadata[order(metadata$condition), ]
metadata$colors <- c("dodgerblue3", "maroon2")[ifelse(metadata$condition == "CTL", 1, 2)]




txdb <- makeTxDbFromGFF(path_gtf, format="gtf")
genes <- genes(txdb)

gtf <- import(path_gtf)

gtf_dexseq <- import(path_gtf_dexseq)
gtf_dexseq <- gtf_dexseq[mcols(gtf_dexseq)$type == "exonic_part"]


mcols(gtf_dexseq)$gene <- mcols(gtf_dexseq)$gene_id
mcols(gtf_dexseq)$id <- mcols(gtf_dexseq)$exonic_part_number


gtf_dexseq


for(j in 1:nrow(valid)){
  # j = 17
  print(j)
  
  g <- valid$gene_id[j]
  chromosome <- as.character(seqnames(genes[g, ]))
  astart <- start(genes[g, ]) - 100
  aend <- end(genes[g, ]) + 100
  
  alTrack <- list()
  sizes <- list()
  
  for(i in 1:nrow(metadata)){
    
    alTrack[[i]] <- AlignmentsTrack(paste0(bam_dir, metadata$sampleName[i], "/accepted_hits.bam"), isPaired = ifelse(metadata$LibraryLayout[i] == "PAIRED", TRUE, FALSE), name = metadata$shortname[i], col = metadata$colors[i], fill = metadata$colors[i], type = c("coverage"))
    sizes[[i]] <- 1
  }
  
  gat <- GenomeAxisTrack()
  alTrack[["gat"]] <- gat
  sizes[["gat"]] <- 1
  
  gtf_sub <- gtf[mcols(gtf)$gene_id == g,]
  
  nr_transcripts <- length(unique(mcols(gtf_sub)$transcript_id))
  
  txTr <- GeneRegionTrack(txdb, name = "transcripts",  transcriptAnnotation = "transcript", just.group = "above", cex.group = 0.7)
  alTrack[["txTr"]] <- txTr
  sizes[["txTr"]] <- ceiling(nr_transcripts/3)
  
  # txTr_dexseq <- GeneRegionTrack(gtf_dexseq[mcols(gtf_dexseq)$gene == g,], name = "htseq bins", exonAnnotation = "exon", collapse = FALSE, fontcolor.exon = 1, cex.exon = 0.7)
  
  gtf_dexseq_sub <- gtf_dexseq[mcols(gtf_dexseq)$gene == g,]
  
  strand(gtf_dexseq_sub) <- "*"
  annTr_dexseq <- AnnotationTrack(gtf_dexseq_sub, name = "htseq bins", showFeatureId = TRUE, fontcolor.feature = "darkblue", stacking = "dense", cex.feature = 0.5)
  
  alTrack[["annTr_dexseq"]] <- annTr_dexseq
  sizes[["annTr_dexseq"]] <- 1
  
  sizes <- unlist(sizes)
  
  pdf(paste0(comparison_out, "gviz/expression_", valid$brooks_gene_id[j], ".pdf"), 12, 8)
  
  try(plotTracks(alTrack, from = astart, to = aend, chromosome = chromosome, sizes = sizes), silent = TRUE)
  
  dev.off()
  
  
  
}












sessionInfo()

