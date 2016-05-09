######################################################
## <<geuvadis_sqtlseeker_2_1.R>>

# BioC 3.1
# Created 25 Nov 2015
# Modified 29 Apr 2016

##############################################################################

Sys.time()

##############################################################################

library(sQTLseekeR)
library(tools)
library(BiocParallel)
library(ggplot2)


##############################################################################
# Arguments for testing the code
##############################################################################

# rwd='/home/Shared/data/seq/geuvadis'
# workers=4
# population='YRI'

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

data_dir <- "data/"

out_dir <- "sqtlseeker_2_1_analysis/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_data_dir <- paste0(out_dir, "data/")
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

out_res_dir <- paste0(out_dir, "results/")
dir.create(out_res_dir, showWarnings = FALSE, recursive = TRUE)

out_plots_dir <- paste0(out_dir, "plots/")
dir.create(out_plots_dir, showWarnings = FALSE, recursive = TRUE)


## Input files: transcript expression, gene location and genotype information
counts_path = paste0(data_dir, paste0("expression/trExpCount_", population, ".tsv"))
gene_bed_path = paste0(data_dir, "annotation/gencode.v12.annotation_genes.bed")
genotypes_path = paste0(out_data_dir, paste0("snps_", population, "_full.tsv"))


##################################################################################
### Prepare genotypes data for sQTLSeekeR
##################################################################################

if(!file.exists(genotypes_path)){
  
  snps_files <- list.files(path = paste0(data_dir, "/genotypes"), pattern = paste0("snps_", population), full.names = TRUE, include.dirs = FALSE)
  
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
  
  cmd <- paste0("cat ", paste0(paste0(out_data_dir, basename(file_path_sans_ext(snps_files)), "_sort.tsv"), collapse = " "), " > ", out_data_dir ,"/snps_", population, ".tsv")
  
  cat(cmd, fill = TRUE)
  
  system(cmd)
  
  
  
  ### add header with samples names
  
  cmd <- paste0("head -n 1 ", snps_files[1], " > ", out_data_dir ,"snps_", population, "_head.tsv")
  
  cat(cmd, fill = TRUE)
  
  system(cmd)
  
  
  cmd <- paste0("cat ", out_data_dir, "snps_", population, "_head.tsv ", out_data_dir, "snps_", population, ".tsv > ", out_data_dir, "snps_", population, "_full.tsv")
  
  cat(cmd, fill = TRUE)
  
  system(cmd)
  
}


###############################################################################
### Run sQTLseekeR analysis 
###############################################################################


## Getting the IDs of samples in the population
metadata <- read.table(paste0(data_dir, "metadata/sample-groups.tsv"), header=TRUE, as.is=TRUE)
metadata <- subset(metadata, group == population)


### check if the samples order is the same in genotype files
stopifnot(all(metadata$sampleShort == read.table(genotypes_path, header = FALSE, nrows = 1)[-c(1:4)]))


########################################
### 1) Index the genotype file (if not done externally before)
########################################

genotypes_path_index <- paste0(out_data_dir, "snps_", population, "_full.tsv.bgz")

if(!file.exists(genotypes_path_index)){
  genotypes_path_index <- index.genotype(genotypes_path)
}


########################################
### 2) Prepare transcript expression
########################################

counts_raw <- read.table(counts_path, as.is=TRUE, header=TRUE, sep="\t")

counts <- counts_raw[,c("trId", "geneId", metadata$sample)]

colnames(counts) <- c("trId", "geneId", metadata$sampleShort)


ratios <- prepare.trans.exp(te.df = counts, min.transcript.exp = 10, min.gene.exp = 10, min.dispersion = 0.01, verbose = FALSE)


write.table(ratios, paste0(out_data_dir, "trExpCount_", population, "_sqtlseeker_ratios.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


ratios <- read.table(paste0(out_data_dir, "trExpCount_", population, "_sqtlseeker_ratios.tsv"), header = TRUE, as.is=TRUE)


tt <- as.numeric(table(ratios$geneId))

df <- data.frame(tt = tt)
binwidth <- ceiling(max(df$tt)/50)

ggp <- ggplot(df, aes_string(x = "tt")) +
  theme_bw() +
  xlab("Number of features per gene") +
  ylab("Frequency") +
  geom_histogram(fill = "seagreen4", binwidth = binwidth) +
  theme(axis.text = element_text(size=16), axis.title = element_text(size=18, face="bold"), plot.title = element_text(size=18, face="bold")) +
  coord_cartesian(xlim = c(0, max(tt) + 2)) +
  geom_text(data = data.frame(x = Inf, y = Inf, label = paste0(length(tt), " genes   \n ", sum(tt) , " features   ")), aes_string(x = "x", y = "y", label = "label"), hjust = 1, vjust = 2, size = 6)


pdf(paste0(out_data_dir, "hist_features_", population, ".pdf"))
print(ggp)
dev.off()


colSums(is.na(counts))
colSums(is.na(ratios))



########################################
### 3) Test gene/SNP associations 
# per block of genes to make it parallel
########################################

# out_res_dir_gene <- paste0(out_res_dir, population, "_gene_level_results/")
# dir.create(out_res_dir_gene)
# 
# 
# gene_bed <- read.table(gene_bed_path, as.is=TRUE, sep="\t")
# colnames(gene_bed) <- c("chr","start","end","geneId")
# 
# 
# genes_unique <- unique(ratios$geneId)
# 
# genes_split <- split(genes_unique, ceiling(seq_along(genes_unique)/100))
# 
# 
# results_list <- bplapply(1:length(genes_split), function(g){
#   # g = 1
#   print(g)
#   
#   res <- sqtl.seeker(ratios[ratios$geneId %in% genes_split[[g]],  , drop = FALSE], genotypes_path_index, gene_bed)
#   
#   if(!is.null(res))
#     write.table(res, paste0(out_res_dir_gene, "results_", g, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#   
#   return(NULL)
#   
# }, BPPARAM = MulticoreParam(workers = workers))
# 
# 
# 
# 
# ##### read in and merge the results 
# 
# results_files <- list.files(out_res_dir_gene, full.names=TRUE)
# 
# results_list <- bplapply(1:length(results_files), function(g){
#   
#   res <- read.table(results_files[g], header = TRUE, as.is = TRUE)
#   
#   return(res)
#   
# }, BPPARAM = MulticoreParam(workers = workers))
# 
# 
# results <- do.call(rbind, results_list)
# 
# 
# write.table(results, paste0(out_res_dir, population, "_results_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

########################################
### 3) Test gene/SNP associations 
# All at once
########################################


gene_bed <- read.table(gene_bed_path, as.is=TRUE, sep="\t")
colnames(gene_bed) <- c("chr","start","end","geneId")

gene_bed$chr <- gsub("chr", "", gene_bed$chr)


results <- sqtl.seeker(tre.df = ratios, genotype.f = genotypes_path_index, gene.loc = gene_bed, genic.window = 5000, min.nb.ext.scores = 1000, nb.perm.max = 1e+06, nb.perm.max.svQTL = 10000, svQTL = FALSE, approx = TRUE, verbose = TRUE)


write.table(results, paste0(out_res_dir, population, "_results_all.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



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



pdf(paste0(out_res_dir, population, "_results_hist_pvalues.pdf"))
print(ggp)
dev.off()




########################################
## 4) Get significant sQTLs
########################################

# results <- read.table(paste0(out_res_dir, population, "_results_all.txt"), header = TRUE, as.is=TRUE)
# 
# 
# results_sign <- sqtls(res.df = results, FDR = 0.05, md.min=.01, out.pdf = paste0(out_res_dir, population, "_results_fdr05.pdf"), svQTL.removal=TRUE, FDR.svQTL=.01)
# 
# 
# write.table(results_sign, paste0(out_res_dir, population, "_results_fdr05.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)









































