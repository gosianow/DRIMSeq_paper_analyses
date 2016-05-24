######################################################
## <<prepare_validated_sqtls.R>>

# BioC 3.2
# Created 1 Mar 2016 
# Updated 9 May 2016

# Process files with validated sQTLs or results from other methods

##############################################################################

Sys.time()

##############################################################################


library(GenomicFeatures)
library(rtracklayer)
library(GenomicRanges)
library(limma)
library(plyr)

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/Shared/data/seq/geuvadis'


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


##########################################################################
### Load data with snp ids and gtf file
##########################################################################

load("data/metadata/snp_id_convert.Rdata")


gtf_file <- "geuvadis_annotation/gencode.v12.annotation.gtf"

gtf0 <- import(gtf_file)

keep <- mcols(gtf0)$type == "gene"
# keep <- mcols(gtf0)$type == "gene" & mcols(gtf0)$gene_type == "protein_coding"

gtf <- gtf0[keep, ]

### Gene ID without dot number
mcols(gtf)$gene_id2 <- strsplit2(mcols(gtf)$gene_id, ".", fixed = TRUE)[, 1]


##########################################################################
### sQTLs from GLIMMPS validated with PCR
##########################################################################

valid_pcr <- read.table("data/validation/glimmps/13059_2013_3131_MOESM4_ESM_processed.txt", sep = "\t", header = TRUE, as.is = TRUE)

### No duplicates
table(duplicated(valid_pcr[, c("SNP.ID", "Ensemble.Gene.ID")]))


match_snps <- match(valid_pcr$SNP.ID, snp_id_convert[, 1])

valid_pcr$snp_id <- snp_id_convert[match_snps, 2]

valid_pcr$snp_name <- valid_pcr$SNP.ID


match_gene <- match(valid_pcr$Gene.symbol, mcols(gtf)$gene_name)

valid_pcr$gene_id <- mcols(gtf)[match_gene, "gene_id"]

valid_pcr$gene_name <- valid_pcr$Gene.symbol


valid_pcr$gene_snp <- paste0(valid_pcr$gene_id, ":", valid_pcr$snp_id)
valid_pcr$chr <- strsplit2(valid_pcr$snp_id, "_")[, 2]
valid_pcr$snp_position <- as.numeric(strsplit2(valid_pcr$snp_id, "_")[, 3])


### Prepare position of PCR target exons
target_exon <- valid_pcr[, "Target.Exon.position..hg19."]

target_exon <- strsplit2(target_exon, ":", fixed = TRUE)[, 2]

target_exon <- strsplit2(target_exon, "-", fixed = TRUE)


valid_pcr$target_exon_start <- target_exon[, 1]
valid_pcr$target_exon_end <- target_exon[, 2]




write.table(valid_pcr, file = "data/validation/glimmps/glimmps_valid_pcr.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




##########################################################################
### sQTLs from GLIMMPS that are related to GWAS SNPs
### Here, we keep the GWAS SNPs as validated associations
##########################################################################

valid_gwas <- read.table("data/validation/glimmps/Table_1_The_list_of_sQTL_signals_linked_to_GWAS_signals_processed.txt", sep = "\t", header = TRUE, as.is = TRUE)

gwas_snp_id <- strsplit2(valid_gwas$GWAS.trait..SNP., split = "(rs", fixed = TRUE)[, 2]

gwas_snp_id <- gsub(")", "", gwas_snp_id, fixed = TRUE)

gwas_snp_id <- paste0("rs", gwas_snp_id)

valid_gwas$gwas_snp_id <- gwas_snp_id

### No duplicates
table(duplicated(valid_gwas[, c("gwas_snp_id", "Gene")]))


match_snps <- match(valid_gwas$gwas_snp_id, snp_id_convert[, 1])

valid_gwas$snp_id <- snp_id_convert[match_snps, 2]

valid_gwas$snp_name <- valid_gwas$gwas_snp_id


match_gene <- match(valid_gwas$Gene, mcols(gtf)$gene_name)

valid_gwas$gene_id <- mcols(gtf)[match_gene, "gene_id"]

valid_gwas$gene_name <- valid_gwas$Gene


valid_gwas$gene_snp <- paste0(valid_gwas$gene_id, ":", valid_gwas$snp_id)
valid_gwas$chr <- strsplit2(valid_gwas$snp_id, "_")[, 2]
valid_gwas$snp_position <- as.numeric(strsplit2(valid_gwas$snp_id, "_")[, 3])


write.table(valid_gwas, file = "data/validation/glimmps/glimmps_valid_gwas.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

##########################################################################
### sQTLs from GLIMMPS that are related to GWAS SNPs
### Here, we keep the GLIMMPS sQTL SNPs as validated associations
##########################################################################

valid_gwas <- read.table("data/validation/glimmps/Table_1_The_list_of_sQTL_signals_linked_to_GWAS_signals_processed.txt", sep = "\t", header = TRUE, as.is = TRUE)

valid_gwas <- unique(valid_gwas[, c("Gene", "sQTL.SNP")])



match_snps <- match(valid_gwas$sQTL.SNP, snp_id_convert[, 1])

valid_gwas$snp_id <- snp_id_convert[match_snps, 2]

valid_gwas$snp_name <- valid_gwas$sQTL.SNP


match_gene <- match(valid_gwas$Gene, mcols(gtf)$gene_name)

valid_gwas$gene_id <- mcols(gtf)[match_gene, "gene_id"]

valid_gwas$gene_name <- valid_gwas$Gene


valid_gwas$gene_snp <- paste0(valid_gwas$gene_id, ":", valid_gwas$snp_id)
valid_gwas$chr <- strsplit2(valid_gwas$snp_id, "_")[, 2]
valid_gwas$snp_position <- as.numeric(strsplit2(valid_gwas$snp_id, "_")[, 3])


write.table(valid_gwas, file = "data/validation/glimmps/glimmps_valid_gwas_glimmps.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



##########################################################################
### SQTLSEEKER 2 validated sQTLs, which in fact come from GLIMMPS
##########################################################################
valid <- read.table("data/validation/sqtlseeker/positive_controls_sqtlseeker_paper.txt", sep = "\t", header = TRUE, as.is = TRUE)


valid$gene_name <- valid$sqtlseeker_gene_id
valid$snp_name <- valid$sqtlseeker_snp_id

valid$gene_snp <- paste0(valid$gene_id, ":", valid$snp_id)
valid$chr <- strsplit2(valid$snp_id, "_")[, 2]
valid$snp_position <- as.numeric(strsplit2(valid$snp_id, "_")[, 3])


write.table(valid, file = "data/validation/sqtlseeker/sqtlseeker_valid.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##########################################################################
### sQTLs found by GLIMMPS 
##########################################################################

valid <- read.table("data/validation/glimmps/13059_2013_3131_MOESM2_ESM_processed.txt", sep = "\t", header = TRUE, as.is = TRUE)


### There are gene-snp pairs that have multiple transcripts significant
table(duplicated(valid[, c("SNP.ID", "Ensembl.Gene.ID")]))

table(duplicated(valid[, c("SNP.ID", "Ensembl.Gene.ID", "Target.Exon.position..hg19.")]))

table(duplicated(valid[, c("SNP.ID", "Ensembl.Gene.ID", "Target.Exon.position..hg19.", "AS.type")]))


valid$snp_name <- valid$SNP.ID

match_snps <- match(valid$snp_name, snp_id_convert[, 1])

valid$snp_id <- snp_id_convert[match_snps, 2]


valid$gene_name <- valid$Gene.Symbol

match_gene <- match(valid$Ensembl.Gene.ID, mcols(gtf)$gene_id2)

valid$gene_id <- mcols(gtf)[match_gene, "gene_id"]


valid$gene_snp <- paste0(valid$gene_id, ":", valid$snp_id)
valid$chr <- strsplit2(valid$snp_id, "_")[, 2]
valid$snp_position <- as.numeric(strsplit2(valid$snp_id, "_")[, 3])


### Prepare positions of target exons
target_exon <- valid[, "Target.Exon.position..hg19."]

target_exon <- strsplit2(target_exon, ":", fixed = TRUE)[, 2]

target_exon <- strsplit2(target_exon, "-", fixed = TRUE)


valid$target_start <- target_exon[, 1]
valid$target_end <- target_exon[, 2]


### Prepare positions of target exons

valid_split <- split.data.frame(valid, factor(valid$gene_snp))


valid_list <- lapply(1:length(valid_split), function(i){
  # i = 4
  valid_tmp <- valid_split[[i]]
  
  if(nrow(valid_tmp) == 1)
    return(valid_tmp)
  
  target_start <- paste0(unique(valid_tmp$target_start), collapse = ",")
  target_end <- paste0(unique(valid_tmp$target_end), collapse = ",")
  
  valid_tmp <- valid_tmp[which.min(valid_tmp$P.value), , drop = FALSE]
  
  valid_tmp$target_start <- target_start
  valid_tmp$target_end <- target_end
  
  return(valid_tmp)
  
})


valid <- rbind.fill(valid_list)


write.table(valid, file = "data/validation/glimmps/glimmps_valid_glimmps.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




##########################################################################
### trQTLs found in GEUVADIS project - all
##########################################################################

valid_set <- data.frame(path = c("EUR373.trratio.cis.FDR5.all.rs137.txt", "YRI89.trratio.cis.FDR5.all.rs137.txt"), population = c("EUR", "YRI"), stringsAsFactors = FALSE)


for(i in 1:nrow(valid_set)){
  
  valid <- read.table(paste0("data/validation/geuvadis/", valid_set$path[i]), sep = "\t", header = TRUE, as.is = TRUE)
  
  
  ### There are gene-snp pairs that have multiple transcripts significant
  table(duplicated(valid[, c("SNP_ID", "GENE_ID")]))
  table(duplicated(valid[, c("SNP_ID", "GENE_ID", "PROBE_ID")]))
  
  table(grepl("snp_", valid$SNP_ID))
  
  
  valid$snp_name <- valid$SNP_ID
  
  match_snps <- match(valid$snp_name, snp_id_convert[, 1])
  
  table(is.na(match_snps))
  
  valid$snp_id <- snp_id_convert[match_snps, 2]
  
  
  valid$gene_id <- valid$GENE_ID
  
  match_gene <- match(valid$gene_id, mcols(gtf)$gene_id)
  
  table(is.na(match_gene))
  
  valid$gene_name <- mcols(gtf)[match_gene, "gene_name"]
  
  
  valid$gene_snp <- paste0(valid$gene_id, ":", valid$snp_id)
  valid$chr <- strsplit2(valid$snp_id, "_")[, 2]
  valid$snp_position <- as.numeric(strsplit2(valid$snp_id, "_")[, 3])
  
  table(duplicated(valid$gene_snp))
  
  
  
  ### There are not only snps in the list but also indels; We keep only snps
  
  valid <- valid[grepl("snp_", valid$snp_id), ]
  
  
  ### Prepare positions of target transcripts and keep only the unique genes with min p-value
  
  valid$target_start <- valid$TSSpos
  valid$target_end <- valid$TSSpos + 1
  
  
  valid_split <- split.data.frame(valid, factor(valid$gene_snp))
  
  
  valid_list <- lapply(1:length(valid_split), function(i){
    # i = 5
    valid_tmp <- valid_split[[i]]
    
    if(nrow(valid_tmp) == 1)
      return(valid_tmp)
    
    PROBE_ID <- paste0(valid_tmp$PROBE_ID, collapse = ",")
    target_start <- paste0(unique(valid_tmp$target_start), collapse = ",")
    target_end <- paste0(unique(valid_tmp$target_end), collapse = ",")
    
    valid_tmp <- valid_tmp[which.min(valid_tmp$pvalue), , drop = FALSE]
    
    valid_tmp$PROBE_ID <- PROBE_ID
    valid_tmp$target_start <- target_start
    valid_tmp$target_end <- target_end
    
    return(valid_tmp)
    
  })
  
  
  valid_geuv <- rbind.fill(valid_list)
  
  
  write.table(valid_geuv, file = paste0("data/validation/geuvadis/geuvadis_valid_geuvadis_all_", valid_set$population[i],".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
}





##########################################################################
### trQTLs found in GEUVADIS project - best
##########################################################################

valid_set <- data.frame(path = c("EUR373.trratio.cis.FDR5.best.rs137.txt", "YRI89.trratio.cis.FDR5.best.rs137.txt"), population = c("EUR", "YRI"), stringsAsFactors = FALSE)


header <- colnames(read.table(paste0("data/validation/geuvadis/YRI89.trratio.cis.FDR5.all.rs137.txt"), sep = "\t", header = TRUE, as.is = TRUE))


for(i in 1:nrow(valid_set)){
  
  valid <- read.table(paste0("data/validation/geuvadis/", valid_set$path[i]), sep = "\t", header = FALSE, as.is = TRUE)
  
  colnames(valid) <- header
  
  
  ### There are gene-snp pairs that have multiple transcripts significant
  table(duplicated(valid[, c("SNP_ID", "GENE_ID")]))
  table(duplicated(valid[, c("SNP_ID", "GENE_ID", "PROBE_ID")]))
  
  table(grepl("snp_", valid$SNP_ID))
  
  
  valid$snp_name <- valid$SNP_ID
  
  match_snps <- match(valid$snp_name, snp_id_convert[, 1])
  
  table(is.na(match_snps))
  
  valid$snp_id <- snp_id_convert[match_snps, 2]
  
  
  valid$gene_id <- valid$GENE_ID
  
  match_gene <- match(valid$gene_id, mcols(gtf)$gene_id)
  
  table(is.na(match_gene))
  
  valid$gene_name <- mcols(gtf)[match_gene, "gene_name"]
  
  
  valid$gene_snp <- paste0(valid$gene_id, ":", valid$snp_id)
  valid$chr <- strsplit2(valid$snp_id, "_")[, 2]
  valid$snp_position <- as.numeric(strsplit2(valid$snp_id, "_")[, 3])
  
  table(duplicated(valid$gene_snp))
  
  
  
  ### There are not only snps in the list but also indels; We keep only snps
  
  valid <- valid[grepl("snp_", valid$snp_id), ]
  
  
  ### Prepare positions of target transcripts and keep only the unique genes with min p-value
  
  valid$target_start <- valid$TSSpos
  valid$target_end <- valid$TSSpos + 1
  
  
  valid_split <- split.data.frame(valid, factor(valid$gene_snp))
  
  
  valid_list <- lapply(1:length(valid_split), function(i){
    # i = 5
    valid_tmp <- valid_split[[i]]
    
    if(nrow(valid_tmp) == 1)
      return(valid_tmp)
    
    PROBE_ID <- paste0(valid_tmp$PROBE_ID, collapse = ",")
    target_start <- paste0(unique(valid_tmp$target_start), collapse = ",")
    target_end <- paste0(unique(valid_tmp$target_end), collapse = ",")
    
    valid_tmp <- valid_tmp[which.min(valid_tmp$pvalue), , drop = FALSE]
    
    valid_tmp$PROBE_ID <- PROBE_ID
    valid_tmp$target_start <- target_start
    valid_tmp$target_end <- target_end
    
    return(valid_tmp)
    
  })
  
  
  valid_geuv <- rbind.fill(valid_list)
  
  
  write.table(valid_geuv, file = paste0("data/validation/geuvadis/geuvadis_valid_geuvadis_best_", valid_set$population[i],".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
}


















