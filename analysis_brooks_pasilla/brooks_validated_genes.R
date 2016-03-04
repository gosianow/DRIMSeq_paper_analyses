######################################################
## ----- brooks_validated_genes
## <<brooks_validated_genes.R>>

# BioC 3.2
# Created 1 Mar 2015 

##############################################################################

Sys.time()

##############################################################################

library(rtracklayer)


##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/brooks_pasilla'
gtf_path='/home/Shared/data/annotation/Drosophila/Ensembl70/gtf/Drosophila_melanogaster.BDGP5.70.gtf'

##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(gtf_path)


##############################################################################

setwd(rwd)


##############################################################################
# validated genes
##############################################################################


valid <- read.table("5_validation/brooks_validated_genes.txt", header = TRUE, sep = "\t", as.is = TRUE) 


##############################################################################
# get gene names 
##############################################################################

gtf <- import(gtf_path)

stopifnot(all(valid$brooks_gene_id %in% mcols(gtf)$gene_name))

keep <- !duplicated(mcols(gtf)$gene_id)
annot <- mcols(gtf)[keep, c("gene_id", "gene_name")]

annot_valid <- annot[annot$gene_name %in% valid$brooks_gene_id, ]

valid$gene_id <- annot_valid[match(valid$brooks_gene_id, annot_valid$gene_name), "gene_id"]
valid$gene_name <- valid$brooks_gene_id


write.table(valid, file = "5_validation/brooks_validated_genes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)






