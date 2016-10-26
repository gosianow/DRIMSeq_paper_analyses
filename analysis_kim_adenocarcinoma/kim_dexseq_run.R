######################################################
## ----- kim_dexseq_run
## <<kim_dexseq_run.R>>

# BioC 2.14
# Created 11 Nov 2015 
# Run DEXSeq_1.10.8

##############################################################################

library("BiocParallel")
library("DEXSeq")

##############################################################################
# Test arguments
##############################################################################

rwd='/home/Shared/data/seq/kim_adenocarcinoma/'
workers=5
count_method=c('htseq','kallisto')[1]
model='model_full2'

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

metadata_org

###############################################################################
# DEXSeq models
###############################################################################


BPPARAM = MulticoreParam(workers = workers)


switch(model,
       
       model_full = {
         
         metadata <- metadata_org
         countFiles <- paste0(counts_out, metadata$sampleName, ".counts")
         sampleTable = data.frame(row.names = metadata$sampleName, condition = metadata$condition)
         formulaFullModel = ~ sample + exon + condition:exon
         formulaReducedModel =  ~ sample + exon 
       },
       
       model_full2 = {
         
         metadata <- metadata_org
         countFiles <- paste0(counts_out, metadata$sampleName, ".counts")
         sampleTable = data.frame(row.names = metadata$sampleName, condition = metadata$condition, Patient.ID = metadata$Patient.ID)
         formulaFullModel = ~ sample + exon + Patient.ID:exon + condition:exon
         formulaReducedModel =  ~ sample + exon + Patient.ID:exon
         
       },
       
       model_null_normal1 = {
         
         metadata <-  metadata_org[metadata_org$condition == "normal", ]
         countFiles <- paste0(counts_out, metadata$sampleName, ".counts")
         sampleTable = data.frame(row.names = metadata$sampleName, condition = rep(c("C1", "C2"), each = 3))
         formulaFullModel = ~ sample + exon + condition:exon
         formulaReducedModel =  ~ sample + exon 
       }, 
       
       model_null_tumor1 = {
         
         metadata <-  metadata_org[metadata_org$condition == "tumor", ]
         countFiles <- paste0(counts_out, metadata$sampleName, ".counts")
         sampleTable = data.frame(row.names = metadata$sampleName, condition = rep(c("C1", "C2"), each = 3))
         formulaFullModel = ~ sample + exon + condition:exon
         formulaReducedModel =  ~ sample + exon 
       },
       
       model_null_normal2 = {
         
         metadata <-  metadata_org[metadata_org$condition == "normal", ]
         countFiles <- paste0(counts_out, metadata$sampleName, ".counts")
         sampleTable = data.frame(row.names = metadata$sampleName, condition = rep(c("C1", "C2"), 3))
         formulaFullModel = ~ sample + exon + condition:exon
         formulaReducedModel =  ~ sample + exon 
       }, 
       
       model_null_tumor2 = {
         
         metadata <-  metadata_org[metadata_org$condition == "tumor", ]
         countFiles <- paste0(counts_out, metadata$sampleName, ".counts")
         sampleTable = data.frame(row.names = metadata$sampleName, condition = rep(c("C1", "C2"), 3))
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



# When using kallisto counts

# Warning messages:
#   1: glm.fit: algorithm did not converge
# 2: glm.fit: algorithm did not converge
# 3: In estimateDispersionsFit(object, fitType = fitType, quiet = quiet) :
#   the parametric fit of dispersion estimates over the mean of counts
# failed, which occurs when the trend is not well captured by the
# function y = a/x + b. A local regression fit is automatically performed,
# and the analysis can continue. You can specify fitType='local' or 'mean'
# to avoid this message if re-running the same data.
# When using local regression fit, the user should examine plotDispEsts(dds)
# to make sure the fitted line is not sharply curving up or down based on
# the position of individual points.













































