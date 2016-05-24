######################################################
## <<brooks_mds.R>>

# BioC 3.2
# Created 2 May 2016 

##############################################################################
Sys.time()
##############################################################################

library(limma)
library(plyr)
library(edgeR)


##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/brooks_pasilla/'
count_method='kallisto'
dir_out='drimseq_0_3_3_comparison'

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

dir_out <- paste0(dir_out, "/")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)


##########################################################################
# load metadata
##########################################################################

metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors = FALSE, sep="\t", header=TRUE) 

metadata_org <- metadata

metadata$colors <- ifelse(metadata$condition == "CTL", "blue", "orange")

metadata$shortname <- paste0(1:nrow(metadata), "_", metadata$condition, "_", substr(metadata$LibraryLayout, 1, 2))

##########################################################################
# Load counts
##########################################################################


count_dir <- paste0("2_counts/", count_method, "/")


### load counts
counts_list <- lapply(1:length(metadata_org$sampleName), function(i){
  # i = 1
  cts <- read.table(paste0(count_dir, metadata_org$sampleName[i], ".txt"), header = FALSE, as.is = TRUE)
  colnames(cts) <- c("group_id", metadata_org$sampleName[i])  
  return(cts)
})

counts <- Reduce(function(...) merge(..., by = "group_id", all=TRUE, sort = FALSE), counts_list)
counts <- counts[!grepl(pattern = "_", counts$group_id),]


group_split <- strsplit2(counts[,1], ":")
counts <- counts[, -1]
### order the samples like in metadata!!!
counts <- counts[, metadata_org$sampleName]




### Calculate gene expression
gene_expression <- by(counts, factor(group_split[, 1]), function(x){
  colSums(x, na.rm = TRUE)
}, simplify = FALSE)

gene_expression <- do.call(rbind, gene_expression)





dge <- DGEList(counts = gene_expression)
dge <- calcNormFactors(dge)
isexpr <- rowSums(cpm(dge) > 1) >= 3
dge <- dge[isexpr, ,keep.lib.sizes=FALSE]



pdf(paste0(dir_out, "pasilla_MDS.pdf"))

mds <- plotMDS(dge, top = 1000, col = metadata$colors, labels = metadata$shortname, cex = 1.5, method="bcv", cex.axis = 1.5, cex.lab = 1.5)

dev.off()

































