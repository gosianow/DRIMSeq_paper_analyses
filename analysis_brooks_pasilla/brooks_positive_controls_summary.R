######################################################
## ----- brooks_positive_controls_summary
## <<brooks_positive_controls_summary.R>>

# BioC 3.1
# Created 30 Nov 2015 

##############################################################################

library(plyr)
library(ggplot2)
library(rtracklayer)
library(reshape2)
library(Gviz)
library(GenomicFeatures)

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/brooks_pasilla'

##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)

##############################################################################

setwd(rwd)

method_out <- "drimseq_0_3_1"

comparison_out <- "drimseq_0_3_1_positive_controls/"
dir.create(comparison_out, showWarnings = FALSE, recursive = TRUE)


gtf_path='/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf'

keep_methods <- c("dexseq", "drimseq_genewise_grid_common", "drimseq_genewise_grid_none")


##############################################################################
# metadata
##############################################################################


metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors=F, sep="\t", header=T) 

metadata

##############################################################################
# validated genes
##############################################################################

valid <- read.table("5_validation/brooks_validated_genes.txt", header = TRUE, sep = "\t", as.is = TRUE) 


##############################################################################
# get gene names 
##############################################################################

gtf <- import(gtf_path)


all(valid$brooks_gene_id %in% mcols(gtf)$gene_name)

keep <- !duplicated(mcols(gtf)$gene_id)
annot <- mcols(gtf)[keep, c("gene_id", "gene_name")]

annot_valid <- annot[annot$gene_name %in% valid$brooks_gene_id, ]

valid$gene_id <- annot_valid[match(valid$brooks_gene_id, annot_valid$gene_name), "gene_id"]

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



summarym <- melt(summary, id.vars = c("model", "count_method"), variable.name = "ds_method", value.name = "counts")
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


pdf(paste0(comparison_out, "validation_summary.pdf"), 15, 5)
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
  
  summarym <- melt(summary[summary$model == models[i], ], id.vars = c("brooks_gene_id", "gene_id", "model", "count_method"), variable.name = "ds_method", value.name = "counts")
  
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
  
  
  pdf(paste0(comparison_out, models[i], "_validation_summary.pdf"), 13, 10)
  print(ggp)
  dev.off()
  
  
}



##############################################################################
# Plot coverage and annotation for the validated genes
##############################################################################

options(ucscChromosomeNames=FALSE)

bam_dir <- "1_reads/tophat_2.0.9/"

metadata <- metadata[order(metadata$condition), ]
metadata$colors <- c("orange", "blueviolet")[ifelse(metadata$condition == "CTL", 1, 2)]


gat <- GenomeAxisTrack()

txdb <- makeTxDbFromGFF(gtf_path, format="gtf")
genes <- genes(txdb)


gene <- valid$gene_id[2]

chromosome <- as.character(seqnames(genes[gene, ]))
astart <- start(genes[gene, ])
aend <- end(genes[gene, ])


txTr <- GeneRegionTrack(txdb, name = "original")



alTrack <- list()

for(i in 1:nrow(metadata)){
  
  alTrack[[i]] <- AlignmentsTrack(paste0(bam_dir, metadata$sampleName[i], "/accepted_hits.bam"), isPaired = ifelse(metadata$LibraryLayout[i] == "PAIRED", TRUE, FALSE), name = metadata$shortname[i], col = metadata$colors[i], fill = metadata$colors[i])
  
}


alTrack[["gat"]] <- gat
alTrack[["txTr"]] <- txTr



pdf("gviz_test.pdf", 10, 10)

plotTracks(alTrack, from = astart, to = aend, chromosome = chromosome, transcriptAnnotation = "transcript", type = c("coverage"), sizes = rep(1/length(alTrack), length(alTrack)))

dev.off()































