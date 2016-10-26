######################################################
## ----- brooks_dexseq_run
## <<brooks_dexseq_run.R>>

# BioC 2.14
# Created 15 Nov 2015 
# Run DEXSeq_1.10.8

##############################################################################

library("BiocParallel")
library("DEXSeq")

##############################################################################
# Test arguments
##############################################################################

# rwd='/home/Shared/data/seq/brooks_pasilla/'
# workers=4
# count_method=c('htseq','kallisto')[2]
# model=c('model_full','model_full2','model_full_paired','model_null1','model_null2','model_null3')[1]


##############################################################################
# Read in the arguments
##############################################################################

## Read input arguments
args <- (commandArgs(trailingOnly = TRUE))
for (i in 1:length(args)) {
  eval(parse(text = args[[i]]))
}


print(rwd)
print(workers)
print(count_method)
print(model)


##############################################################################

setwd(rwd)

method_out <- "4_results/dexseq_1_10_8/"

counts_out <- paste0("2_counts/", count_method, "/")

out_dir <- paste0(method_out, "/",  model, "/", count_method)
dir.create(out_dir, recursive = TRUE)


###############################################################################
# metadata
###############################################################################


metadata <- read.table("3_metadata/metadata.xls", stringsAsFactors = FALSE, sep="\t", header=TRUE) 

metadata_org <- metadata


###############################################################################
# DEXSeq models
###############################################################################


BPPARAM = MulticoreParam(workers = workers)


switch(model,
       
       model_full = {
         
         metadata <- metadata_org
         countFiles <- paste0(counts_out, metadata$SampleName, ".txt")
         sampleTable = data.frame(row.names = metadata$SampleName, condition = metadata$condition)
         formulaFullModel = ~ sample + exon + condition:exon
         formulaReducedModel =  ~ sample + exon 
         
       },
       
       model_full2 = {
         
         metadata <- metadata_org
         countFiles <- paste0(counts_out, metadata$SampleName, ".txt")
         sampleTable = data.frame(row.names = metadata$SampleName, condition = metadata$condition, LibraryLayout = metadata$LibraryLayout)
         formulaFullModel = ~ sample + exon + LibraryLayout:exon + condition:exon
         formulaReducedModel =  ~ sample + exon + LibraryLayout:exon
         
       }, 
       
       model_full_paired = {
         
         metadata <-  metadata_org[metadata_org$LibraryLayout == "PAIRED", ]
         countFiles <- paste0(counts_out, metadata$SampleName, ".txt")
         sampleTable = data.frame(row.names = metadata$SampleName, condition = metadata$condition)
         formulaFullModel = ~ sample + exon + condition:exon
         formulaReducedModel =  ~ sample + exon
         
       },
       
       model_null1 = {
         
         metadata <-  metadata_org[metadata_org$condition == "CTL", ]
         countFiles <- paste0(counts_out, metadata$SampleName, ".txt")
         sampleTable = data.frame(row.names = metadata$SampleName, condition = rep(c("C1", "C2"), 2))
         formulaFullModel = ~ sample + exon + condition:exon
         formulaReducedModel =  ~ sample + exon
         
       }, 
       
       model_null2 = {
         
         metadata <-  metadata_org[metadata_org$condition == "CTL", ]
         countFiles <- paste0(counts_out, metadata$SampleName, ".txt")
         sampleTable = data.frame(row.names = metadata$SampleName, condition = rep(c("C1", "C2"), each = 2))
         formulaFullModel = ~ sample + exon + condition:exon
         formulaReducedModel =  ~ sample + exon
         
       },
       
       model_null3 = {
         
         metadata <-  metadata_org[metadata_org$condition == "CTL", ]
         countFiles <- paste0(counts_out, metadata$SampleName, ".txt")
         sampleTable = data.frame(row.names = metadata$SampleName, condition = c("C1", "C2", "C2", "C1"))
         formulaFullModel = ~ sample + exon + condition:exon
         formulaReducedModel =  ~ sample + exon
         
       }
       
)


dxd <- DEXSeqDataSetFromHTSeq(countFiles, sampleData = sampleTable, design = ~ sample + exon + condition:exon , flattenedfile = NULL ) 

dxr <- DEXSeq(dxd, fullModel = formulaFullModel, reducedModel = formulaReducedModel, BPPARAM = BPPARAM, fitExpToVar="condition")


dxrf <- as.data.frame(dxr)

write.table(dxrf[, -ncol(dxrf)], paste0(out_dir, "/dexseq_exon_results.txt"), sep="\t", row.names=FALSE, quote=F) 

### gene level test
perGeneRes <- perGeneQValue(dxr)

write.table(data.frame(geneID = names(perGeneRes), pvalue = perGeneRes), paste0(out_dir, "/dexseq_gene_results.txt"), sep="\t", row.names=FALSE, quote=F)





sessionInfo()





