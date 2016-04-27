######################################################
## <<geuvadis_sqtlseeker_2_1_drimseq_counts.R>>

# BioC 3.1
# Created 25 Nov 2015
# Modified 26 Apr 2016

##############################################################################

Sys.time()

##############################################################################

library(sQTLseekeR)
library(tools)
library(BiocParallel)
library(ggplot2)
library(plyr)

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/Shared/data/seq/geuvadis'
workers=10
drimseq_results_path='drimseq_0_3_3_analysis_permutations_all_genes'

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
print(workers)


##############################################################################

setwd(rwd)

data_dir <- "data/"

out_dir <- "sqtlseeker_2_1_analysis_drimseq_counts/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_data_dir <- paste0(out_dir, "data/")
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

out_res_dir <- paste0(out_dir, "results/")
dir.create(out_res_dir, showWarnings = FALSE, recursive = TRUE)

out_plots_dir <- paste0(out_dir, "plots/")
dir.create(out_plots_dir, showWarnings = FALSE, recursive = TRUE)


## Input files: transcript expression, gene location and genotype information
gene_bed_path = paste0(data_dir, "annotation/gencode.v12.annotation_genes.bed")
genotypes_path = paste0(out_data_dir, "snps_CEU_full.tsv")


##################################################################################
### Prepare genotypes data for sQTLSeekeR
##################################################################################

if(!file.exists(genotypes_path)){
  
  snps_files <- list.files(path = paste0(data_dir, "/genotypes"), pattern = "snps_CEU", full.names = TRUE, include.dirs = FALSE)
  
  snps_files <- snps_files[grepl("chr",snps_files) & !grepl("sort",snps_files)]
  
  x <- gsub("^.*chr","", snps_files)
  xx <- gsub(".tsv$", "", x)
  
  snps_files <- snps_files[order(as.numeric(xx))]
  
  snps_files
  
  ### sort by SNP position
  
  for(i in 1:length(snps_files)){
    # i = 1
    cat(i, fill = TRUE)
    #### Does not sort in the right way -->> use tab as separator
    #   cmd <- paste0("(head -n 1 ", snps_files[i], " && tail -n +2 ", snps_files[i], " | sort -k 2 ) > ", snps_files[i], ".sort.tsv")
    #   cmd <- paste0( "tail -n +2 ", snps_files[i], " | sort -k 2  > ", snps_files[i], ".sort.tsv")  
    #   system(cmd)
    
    d <- read.table(snps_files[i], header = TRUE, as.is = TRUE)
    
    o <- order(d[,2])
    
    print(table(o[-1] - o[-length(o)]))
    
    do <- d[o, ]
    
    write.table(do, paste0(out_data_dir, basename(file_path_sans_ext(snps_files[i])), "_sort.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    
    gc()
    
  }
  
  
  ### merge chromosome files
  
  cmd <- paste0("cat ", paste0(paste0(out_data_dir, basename(file_path_sans_ext(snps_files)), "_sort.tsv"), collapse = " "), " > ", out_data_dir ,"/snps_CEU.tsv")
  
  cat(cmd, fill = TRUE)
  
  system(cmd)
  
  
  
  ### add header with samples names
  
  cmd <- paste0("head -n 1 ", snps_files[1], " > ", out_data_dir ,"snps_CEU_head.tsv")
  
  cat(cmd, fill = TRUE)
  
  system(cmd)
  
  
  cmd <- paste0("cat ", out_data_dir, "snps_CEU_head.tsv ", out_data_dir, "snps_CEU.tsv > ", out_data_dir, "snps_CEU_full.tsv")
  
  cat(cmd, fill = TRUE)
  
  system(cmd)
  
  
}



###############################################################################
### Run sQTLseekeR analysis 
###############################################################################




## Getting the IDs of samples in CEU population
metadata <- read.table(paste0(data_dir, "metadata/sample-groups.tsv"), header=TRUE, as.is=TRUE)
metadata <- subset(metadata, group == "CEU")


### check if the samples order is the same in genotype files
stopifnot(metadata$sampleShort == read.table(genotypes_path, header = FALSE, nrows = 1)[-c(1:4)])


########################################
### 1) Index the genotype file (if not done externally before)
########################################

genotypes_path_index <- index.genotype(genotypes_path)
genotypes_path_index
# genotypes_path_index <- paste0(out_data_dir, "snps_CEU_full.tsv.bgz")


########################################
### 2) Prepare transcript expression
########################################

### Use prefiltered DRIMSeq counts for the analysis 


counts <- lapply(1:22, function(chr){
  # chr = 1
  
  load(paste0(drimseq_results_path, "/CEU_chr",chr, "_d.Rdata"))
  
  out <- data.frame(trId = rownames(d@counts), geneId = rep(names(d@counts@partitioning), times = width(d@counts)), d@counts@unlistData, stringsAsFactors = FALSE)
  
  return(out)
  
})

counts <- rbind.fill(counts)

counts <- counts[,c("trId", "geneId", metadata$sampleShort)]


ratios <- by(counts[, metadata$sampleShort], factor(counts$geneId, levels = unique(counts$geneId)), function(x){
  x <- as.matrix(x)
  out <- data.frame(sweep(x, 2, colSums(x, na.rm = TRUE), "/"))
  return(out)
}, simplify = FALSE)


ratios <- rbind.fill(ratios)
ratios <- data.frame(trId = counts$trId, geneId = counts$geneId, ratios, stringsAsFactors = FALSE)


write.table(ratios, paste0(out_data_dir, "trExpCount_CEU_sqtlseeker_ratios.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


ratios <- read.table(paste0(out_data_dir, "trExpCount_CEU_sqtlseeker_ratios.tsv"), header = TRUE, as.is=TRUE)


tt <- as.numeric(table(ratios$geneId))

df <- data.frame(tt = tt)
binwidth <- 1

ggp <- ggplot(df, aes_string(x = "tt")) +
  theme_bw() +
  xlab("Number of features per gene") +
  ylab("Frequency") +
  geom_histogram(fill = "seagreen4", binwidth = binwidth) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold")) +
  coord_cartesian(xlim = c(0, max(tt) + 2)) +
  geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(length(tt), " genes   \n ", sum(tt) , " features   ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, size = 6)


pdf(paste0(out_data_dir, "hist_features.pdf"))
print(ggp)
dev.off()


colSums(is.na(counts))
colSums(is.na(ratios))






########################################
### 3) Test gene/SNP associations 
# per block of genes to make it parallel
########################################

out_res_dir_gene <- paste0(out_res_dir, "CEU_gene_level_results/")
dir.create(out_res_dir_gene)


gene_bed <- read.table(gene_bed_path, as.is=TRUE, sep="\t")
colnames(gene_bed) <- c("chr","start","end","geneId")


genes_unique <- unique(ratios$geneId)

genes_split <- split(genes_unique, ceiling(seq_along(genes_unique)/100))


results_list <- bplapply(1:length(genes_split), function(g){
  # g = 1
  print(g)
  
  res <- sqtl.seeker(ratios[ratios$geneId %in% genes_split[[g]],  , drop = FALSE], genotypes_path_index, gene_bed)
  
  if(!is.null(res))
    write.table(res, paste0(out_res_dir_gene, "results_", g, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  return(NULL)
  
}, BPPARAM = MulticoreParam(workers = workers))




##### read in and merge the results 

results_files <- list.files(out_res_dir_gene, full.names=TRUE)

results_list <- bplapply(1:length(results_files), function(g){
  
  res <- read.table(results_files[g], header = TRUE, as.is = TRUE)
  
  return(res)
  
}, BPPARAM = MulticoreParam(workers = workers))


results <- do.call(rbind, results_list)


write.table(results, paste0(out_res_dir, "CEU_results_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

results <- read.table(paste0(out_res_dir, "CEU_results_all.txt"), header = TRUE, as.is=TRUE)

########################################
## 4) Get significant sQTLs
########################################


results_sign = sqtls(results, FDR = 0.05, out.pdf = paste0(out_res_dir, "CEU_results_fdr05.pdf"))


write.table(results_sign, paste0(out_res_dir, "CEU_results_fdr05.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


pvalues <- results[, "pv"]

df <- data.frame(pvalues = pvalues[!is.na(pvalues)])

ggp <- ggplot(df, aes_string(x = "pvalues")) +
  theme_bw() +
  xlab("p-values") +
  ylab("Frequency") +
  geom_histogram(breaks = seq(0, 1, by = 0.01), fill = "deeppink4") +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=16, face="bold")) +
  coord_cartesian(xlim = c(-0.02, 1.02)) +
  geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(nrow(df), " tests       ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 3, size = 6)



pdf(paste0(out_res_dir, "CEU_results_fdr05_", "hist_features.pdf"))
print(ggp)
dev.off()









































