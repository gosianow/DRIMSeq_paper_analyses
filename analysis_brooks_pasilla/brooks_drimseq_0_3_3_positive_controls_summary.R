######################################################
## ----- brooks_drimseq_0_3_3_positive_controls_summary
## <<brooks_drimseq_0_3_3_positive_controls_summary.R>>

# BioC 3.2
# Created 15 Jan 2015 

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

method_out <- "drimseq_0_3_3"

comparison_out <- "drimseq_0_3_3_positive_controls/"
dir.create(comparison_out, showWarnings = FALSE, recursive = TRUE)

dir.create(paste0(comparison_out, "figures/"), showWarnings = FALSE, recursive = TRUE)

path_gtf_dexseq <- paste0(file_path_sans_ext(path_gtf), ".DEXSeq.flattened.rNO.gff")
path_gtf_filtered_dexseq <- paste0(file_path_sans_ext(path_gtf_filtered), ".DEXSeq.flattened.rNO.gff")


keep_methods <- c("dexseq", "drimseq_genewise_grid_common", "drimseq_genewise_grid_none")


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
# merge positive controls summary into one table
##############################################################################


files <- list.files(comparison_out, pattern = "_validation_summary.txt")
files


summary_list <- lapply(1:length(files), function(i){
  
  sm <- read.table(paste0(comparison_out, files[i]), header = TRUE, as.is = TRUE)
  
})


summary <- rbind.fill(summary_list)

summary <- summary[, c("model", "count_method", keep_methods)]


write.table(summary, file = paste0(comparison_out, "validation_summary.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



summarym <- melt(summary, id.vars = c("model", "count_method"))
# summarym <- melt(summary, id.vars = c("model", "count_method"), variable.name = "ds_method", value.name = "counts")

summarym$variable <- factor(summarym$variable, levels = rev(keep_methods))


ggp <- ggplot(summarym, aes(x = count_method, y = variable, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = value), color = "black", size = 7) + 
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) + 
  xlab("") + 
  ylab("") + 
  theme_bw() +
  theme(panel.background = element_rect(fill = NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14), strip.text = element_text(size = 14)) +
  scale_fill_gradient(low = "grey90", high = "grey50", na.value = "grey90") +
  facet_wrap(~ model)


pdf(paste0(comparison_out, "figures/validation_summary.pdf"), 15, 5)
print(ggp)
dev.off()



##############################################################################
# merge positive controls into tables per model
##############################################################################


files <- list.files(comparison_out, pattern = "_validation.txt")
files


summary_list <- lapply(1:length(files), function(i){
  
  sm <- read.table(paste0(comparison_out, files[i]), header = TRUE, as.is = TRUE)
  
})


summary <- rbind.fill(summary_list)

summary <- summary[, c("brooks_gene_id", "gene_id", "model", "count_method", keep_methods)]
summary$model <- factor(summary$model)

models <- levels(summary$model)




for(i in 1:nlevels(summary$model)){
  # i = 1
  
  # summarym <- melt(summary[summary$model == models[i], ], id.vars = c("brooks_gene_id", "gene_id", "model", "count_method"), variable.name = "ds_method", value.name = "counts")
  summarym <- melt(summary[summary$model == models[i], ], id.vars = c("brooks_gene_id", "gene_id", "model", "count_method"))
  
  summarym$variable <- factor(summarym$variable, levels = keep_methods)
  
  summarym$status <- factor(summarym$value < 0.05)
  summarym$full_gene_id <- paste0(summarym$brooks_gene_id, "_", summarym$gene_id)
  
  ggp <- ggplot(summarym, aes(x = variable, y = full_gene_id, fill = status)) + 
    geom_tile() + 
    geom_text(aes(label = sprintf( "%.02e", value)), color = "black", size = 4) + 
    scale_x_discrete(expand = c(0, 0)) + 
    scale_y_discrete(expand = c(0, 0)) + 
    xlab("") + 
    ylab("") + 
    theme_bw() +
    theme(panel.background = element_rect(fill = NA, colour = NA), axis.ticks = element_blank(), axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(size = 14), strip.text = element_text(size = 14)) +
    scale_fill_manual(values = c("grey80", "grey50"), na.value = "grey90") +
    facet_wrap(~ count_method, nrow = 1)
  
  
  pdf(paste0(comparison_out, "figures/", models[i], "_validation_summary.pdf"), 13, 10)
  print(ggp)
  dev.off()
  
  
}



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




















