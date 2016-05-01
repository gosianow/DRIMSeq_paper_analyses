##############################################################################
### Prepare data 
##############################################################################


# BioC 3.0
# Created 9 Jan 2014
# Updated 1 May 2016


##############################################################################

setwd("/home/Shared/data/seq/geuvadis/")

out_dir <- "data/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

##############################################################################
###  gene.bed.f
##############################################################################

library(GenomicRanges)
library(rtracklayer)

gtfFile <- "geuvadis_annotation/gencode.v12.annotation.gtf"

gtf0 <- import(gtfFile)

## keep protein coding genes
pcindex <- mcols(gtf0)$gene_type == "protein_coding" 
gtf <- gtf0[pcindex]

gindex <- mcols(gtf)$type == "gene" 
gtfg <- gtf[gindex]

gene.bed.f <- data.frame(chr = seqnames(gtfg), start =  start(gtfg), end = end(gtfg), geneId = mcols(gtfg)$gene_id)

write.table(gene.bed.f, paste0(out_dir,"annotation/gencode.v12.annotation_genes.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



# gene.bed <- read.table(paste0(out_dir,"genes.bed"))
# 
# gene.bed[,1] <- gsub(pattern = "chr", replacement = "", x = gene.bed[,1])
# 
# write.table(gene.bed, paste0(out_dir,"genes_noChr.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)



##############################################################################
###  metadata
##############################################################################

metadata <- read.table("geuvadis_analysis_results/E-GEUV-1.sdrf.txt", header = T, sep="\t", as.is=TRUE)

metadata <- metadata[c("Assay.Name", "Characteristics.population.")]
metadata <- unique(metadata)

table(metadata$Characteristics.population.)

colnames(metadata) <- c("sample", "group")
metadata$sampleShort <- substr(metadata$sample, 1, 7)

write.table(metadata, paste0(out_dir, "metadata/sample-groups.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE) 




##############################################################################
###  trans.exp.f
##############################################################################

### RPKM
quantFile <- "geuvadis_analysis_results/GD660.TrQuantRPKM.txt"

trans.exp.f0 <- read.table(quantFile, header = T, sep="\t", as.is = TRUE)

trans.exp.f <- trans.exp.f0[, c("TargetID", "Gene_Symbol", metadata$sample)]
colnames(trans.exp.f) <- c("trId", "geneId", metadata$sample)

write.table(trans.exp.f, paste0(out_dir,"expression/trExpRPKM.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

for(i in unique(metadata$group))
  write.table(trans.exp.f[,c("trId", "geneId", metadata$sample[metadata$group == i])], paste0(out_dir,"expression/trExpRPKM_",i,".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



### Counts
quantFile <- "geuvadis_analysis_results/GD660.TrQuantCount.txt"

trans.exp.f0 <- read.table(quantFile, header = T, sep="\t", as.is = TRUE)

trans.exp.f <- trans.exp.f0[, c("TargetID", "Gene_Symbol", metadata$sample)]
colnames(trans.exp.f) <- c("trId", "geneId", metadata$sample)

write.table(trans.exp.f, paste0(out_dir,"expression/trExpCount.tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

for(i in unique(metadata$group))
  write.table(trans.exp.f[,c("trId", "geneId", metadata$sample[metadata$group == i])], paste0(out_dir,"expression/trExpCount_",i,".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




# ### Check the relation between counts and RPKM
# 
# trans.exp.rpkm <- read.table(paste0(out_dir,"expression/trExpRPKM_CEU.tsv"), header = T, sep="\t", as.is = TRUE)
# 
# trans.exp.counts <- read.table(paste0(out_dir,"expression/trExpCount_CEU.tsv"), header = T, sep="\t", as.is = TRUE)
# 
# all(trans.exp.rpkm$trId == trans.exp.counts$trId)
# 
# 
# pdf(paste0(out_dir,"expression/Relation_between_RPKM&Count.pdf"))
# plot(trans.exp.rpkm[,3], trans.exp.counts[,3])
# dev.off()


##############################################################################
###  genotype.f 
##############################################################################

library(VariantAnnotation)
library(limma)


tbindex <- list.files(path = "geuvadis_genotypes", pattern = "genotypes.vcf.bgz.tbi", full.names = TRUE, include.dirs = FALSE)
compressVcf <- substr(tbindex, 1, nchar(tbindex)-4)


chr <- gsub("chr", "", strsplit2(compressVcf, split=".", fixed=TRUE)[,2])


### extended gene ranges
gnrng <- gtfg 
range <- 5000
start(gnrng) <- start(gnrng) - range
end(gnrng) <- end(gnrng) + range
seqlevels(gnrng) <- gsub("chr", "", seqlevels(gnrng))




for(j in c("CEU", "FIN", "GBR", "TSI", "YRI")){
  
  for(i in 1:length(compressVcf)){
    # j = "TSI"; i = 3
    
    message(paste0(j, chr[i]))
    
    tab <- TabixFile(compressVcf[i], tbindex[i])
    gnrngTmp <- gnrng[seqnames(gnrng) == chr[i]]
    sum(width(gnrngTmp))
    
    gnrngTmp <- reduce(gnrngTmp)
    sum(width(gnrngTmp))
    
    ## Explore the file header with scanVcfHeader
    hdr <- scanVcfHeader(tab)
    
    stopifnot(all(metadata$sampleShort %in% samples(hdr)))

    
    ## read VCF file 
    param <- ScanVcfParam(which = gnrngTmp, samples = metadata$sampleShort[metadata$group == j])
    
    vcf <- readVcf(tab, "hg19", param)
    #   vcf
    
    
    ## Keep only bi-allelic SNPs
    
    # width of ref seq
    rw <- width(ref(vcf))
    # width of first alt seq
    aw <- unlist(lapply(alt(vcf), function(x) {width(x[1])}))
    # number of alternate genotypes
    nalt <- elementLengths(alt(vcf))
    # select only bi-allelic SNPs (monomorphic OK, so aw can be 0 or 1)
    snp <- rw == 1 & aw <= 1 & nalt == 1
    
    # subset vcf  
    vcfbi <- vcf[snp,]
    # remove duplicates
    vcfbi <- vcfbi[!duplicated(rownames(vcfbi)), ]
    
    
    rowdata <- rowData(vcfbi)
    
    ## Convert genotype into number of alternative allels
    geno <- geno(vcfbi)$GT
    
    geno01 <- matrix(-1, nrow = nrow(geno), ncol = ncol(geno))
    rownames(geno01) <- rownames(geno)
    colnames(geno01) <- colnames(geno)
    
    geno01[geno %in% c("0/0", "0|0")] <- 0 # REF/REF
    geno01[geno %in% c("0/1", "0|1", "1/0", "1|0")] <- 1 # REF/ALT
    geno01[geno %in% c("1/1", "1|1")] <- 2 # ALT/ALT
    
    ### genotype
    genotype <- data.frame(chr = seqnames(rowdata), start = start(rowdata), end = end(rowdata), snpId = rownames(geno01), geno01)
    
    
    ## ref alt
    ref <- as.matrix(ref(vcfbi))
    alt <- unlist(alt(vcfbi))
    
    ref_alt <- data.frame(chr = seqnames(rowdata), start = start(rowdata), end = end(rowdata), snpId = rownames(geno01), ref = ref, alt = alt)
    
    
    ### sorting
    oo <- order(genotype[,2])
    
    genotype <- genotype[oo, ]
    ref_alt <- ref_alt[oo, ]
    
    
    write.table(genotype, file=paste0(out_dir, "genotypes/snps_", j, "_chr", chr[i], ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    write.table(ref_alt, file=paste0(out_dir, "genotypes/refalt_", j, "_chr", chr[i], ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    
    
    gc()
    
  }
  
  
}



##############################################################################
### Process SNP ID data
### Import Phase1.Geuvadis_dbSnp137_idconvert_snpOnly.txt
##############################################################################


snp_id_convert <- read.table(paste0("geuvadis_genotypes/Phase1.Geuvadis_dbSnp137_idconvert.txt"), header = FALSE, as.is= TRUE)

save(snp_id_convert, file = paste0(out_dir, "metadata/snp_id_convert.Rdata"))











































