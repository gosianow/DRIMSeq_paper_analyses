
# R32dev

library(devtools)
load_all("/home/gosia/R/multinomial_project/package_devel/DRIMSeq")

##############################################################################
# GEUVADIS 
##############################################################################

library(DRIMSeq)
library(ggplot2)
library(limma)
library(GenomicRanges)
library(plyr)
library(reshape2)

##############################################################################
# Arguments for testing the code
##############################################################################

rwd='/home/Shared/data/seq/geuvadis'
workers=5
population='CEU'
chr='22'


setwd(rwd)

out_dir <- "drimseq_0_3_3_analysis_f/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_name <- paste0(out_dir, population, "_chr",chr, "_")


### Results with CR adjustment
load(paste0(out_name, "d_cr.Rdata"))

dcr <- d

### Results with NO CR adjustment
load(paste0(out_name, "d_ncr.Rdata"))

dncr <- d


### New output directory
out_dir <- "drimseq_0_3_3_analysis_f/investigate_results/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_name <- paste0(out_dir, population, "_chr", chr, "_")


#####################################
### sqtlseeker results
#####################################


resseeker <- read.table(paste0("sqtlseeker_2_1_analysis/results/", population, "_results_all.txt"), header = TRUE, as.is = TRUE)
head(resseeker)

colnames(resseeker) <- c("gene_id", "snp_id", "F", "nb.groups", "md", "tr.first", "tr.second", "nb.perms", "pvalue")
resseeker <- unique(resseeker)


resseeker$gene_snp <- paste0(resseeker$gene_id, ":", resseeker$snp_id)

resseeker$adj_pvalue <- qvalue::qvalue(resseeker$pvalue)$qvalues

### Keep results for chr of interest
resseeker <- resseeker[grepl(paste0("snp_", chr), resseeker$snp_id), ]


#####################################
### merge results
#####################################



rescr <- results(dcr)
resncr <- results(dncr)

rescr$gene_block <- paste0(rescr$gene_id, ":", rescr$block_id)
rescr$gene_snp <- paste0(rescr$gene_id, ":", rescr$snp_id)

resncr$gene_block <- paste0(resncr$gene_id, ":", resncr$block_id)
resncr$gene_snp <- paste0(resncr$gene_id, ":", resncr$snp_id)



table(rescr$adj_pvalue < 0.05)
table(resncr$adj_pvalue < 0.05)
table(resseeker$adj_pvalue < 0.05)


colnames(rescr) <- paste0(colnames(rescr), "_cr")
colnames(resncr) <- paste0(colnames(resncr), "_ncr")
colnames(resseeker) <- paste0(colnames(resseeker), "_seeker")

resm <- merge(rescr, resncr, by.x = "gene_snp_cr", by.y = "gene_snp_ncr", sort = FALSE)

resm <- merge(resm, resseeker, by.x = "gene_snp_cr", by.y = "gene_snp_seeker", sort = FALSE)




#####################################
### Plot some examples and simulate similar scenarios
#####################################

simulation_script='/home/gosia/R/drimseq_paper/simulations_dm/dm_simulate.R'
source(simulation_script)


### Check if sqtlseeker gives the same results for snps from one block - YES
all_equal <- by(resm, resm$gene_block_cr, function(x){
  sd(x$pvalue_seeker, na.rm = TRUE) 
})

table(all_equal, useNA = "always")





### Plot examples where significant for DRIMSEQ and not significant for SQTLSEEKER 

resm_sign <- resm[resm$adj_pvalue_cr < 0.05, ]

resm_sign <- resm_sign[order(resm_sign$pvalue_seeker, decreasing = TRUE), ]

genes <- resm_sign[, c("gene_id_cr", "snp_id_cr", "block_id_cr")]
colnames(genes) <- c("gene_id", "snp_id", "block_id")

genes <- genes[!duplicated(genes[, c("gene_id")]), ]



for(i in 1:10){
  # i = 1
  plotFit(dcr, gene_id = genes[i, "gene_id"], snp_id = genes[i, "snp_id"], out_dir = paste0(out_name, "dm0seeker1_cr_", i, "_"), order = FALSE)
  
  genotype <- dcr@genotypes[[genes[i, "gene_id"]]][genes[i, "block_id"], ]
  genotype <- genotype[!is.na(genotype)]
  pi <- dcr@fit_null[[genes[i, "gene_id"]]][[genes[i, "block_id"]]][, "null"]
  g0 <- dcr@genewise_dispersion[[genes[i, "gene_id"]]][genes[i, "block_id"]]
  nm <- dcr@mean_expression[genes[i, "gene_id"]]
  
  
  counts <- dm_simulate(m = 1, n = length(genotype), pi = pi, g0 = g0, nm = nm, nd = 0, mc.cores = workers)
  
  group_split <- strsplit2(rownames(counts), ":")
  group_split[, 1] <- genes[i, "gene_id"]
  group_split[, 2] <- names(pi)
  
  ### Treat data like for DS analysis
  d <- dmDSdata(counts = counts, gene_id = group_split[, 1], feature_id = group_split[, 2], sample_id = paste0("s", 1:ncol(counts)), group = factor(genotype))
  
  d <- dmDispersion(d, mean_expression = FALSE, common_dispersion = FALSE, genewise_dispersion = TRUE, disp_adjust = FALSE, disp_init = g0)
  
  d <- dmFit(d)
  
  d <- dmTest(d)
  
  plotFit(d, gene_id = genes[i, "gene_id"], out_dir = paste0(out_name, "dm0seeker1_cr_", i, "_simulation_"), plot_type = "boxplot1", order = FALSE)
  
}



for(i in 1:10){
  
  plotFit(dncr, gene_id = genes[i, "gene_id"], snp_id = genes[i, "snp_id"], out_dir = paste0(out_name, "dm0seeker1_ncr_", i, "_"), order = FALSE)
  
}



### Check the contribution of each transcript to the LR test

lik_contribution <- function(y, pi, gamma0){
  
  l <- rowSums(lgamma(y + pi * gamma0) - lgamma(pi * gamma0), na.rm = TRUE)
  
  return(l)
  
}



# d = dcr; gene_id = genes[i, "gene_id"]; snp_id = genes[i, "snp_id"]; name = paste0(out_name, "dm0seeker1_cr_", i, "_")

lr_contribution <- function(d, gene_id, snp_id, name){
  
  block_id <- d@blocks[[gene_id]][d@blocks[[gene_id]][, "snp_id"] == snp_id, "block_id"]
  
  gamma0 <- slot(d, d@dispersion)
  
  if(d@dispersion == "genewise_dispersion")
    gamma0 <- gamma0[[gene_id]][block_id]
  
  genotypes <- d@genotypes[[gene_id]][block_id, ]
  counts <- d@counts[[gene_id]]
  NAs <- !(is.na(counts[1,]) | is.na(genotypes))
  
  ## Null
  y <- counts[, NAs, drop = FALSE]
  pi <- d@fit_null[[gene_id]][[block_id]][, "null"]
  
  likn <- lik_contribution(y, pi, gamma0)
  
  ### Check if calculation is correct
  N <- ncol(y)
  S <- colSums(y)
  l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + sum(likn, na.rm = TRUE)
  
  d@fit_null[[gene_id]]@metadata[block_id, "lik"] - l
  
  
  ### Full
  genotypes <- factor(genotypes[NAs])
  counts <- counts[, NAs, drop = FALSE]
  groups <- levels(genotypes)
  
  likf <- sapply(1:length(groups), function(ii){
    # ii = 2
    
    group <- groups[ii]
    y <- counts[, genotypes == group]
    pi <- d@fit_full[[gene_id]][[block_id]][, group]
    
    likfg <- lik_contribution(y, pi, gamma0)
    
    ### Check if calculation is correct
    N <- ncol(y)
    S <- colSums(y)
    l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + sum(likfg, na.rm = TRUE)
    
    d@fit_full[[gene_id]]@metadata[block_id, paste0("lik_", group)] - l
    
    return(likfg)
    
  })
  
  lrc <- 2 * (rowSums(likf, na.rm = TRUE) - likn)
  
  lrc <- data.frame(t(lrc))
  
  lrc$gene_id <- gene_id
  lrc$block_id <- block_id
  lrc$snp_id <- snp_id
  
  lrcm <- melt(lrc, id.vars = c("gene_id", "block_id", "snp_id"), variable.name = "feature", value.name = "lr_contrib")
  
  lrcm$feature <- factor(lrcm$feature, levels = unique(lrcm$feature))
  
  lrc_tmp <- lrcm
  
  lr_tmp <- sum(lrc_tmp$lr_contrib, na.rm = TRUE)
  
  ggp <- ggplot(lrc_tmp, aes(x = feature, y = lr_contrib)) +
    geom_bar(stat = "identity") +
    ggtitle(paste0("LR = ", lr_tmp)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16), axis.text.y = element_text(size = 16), axis.title.x = element_blank(), axis.title.y = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), legend.title = element_blank(), legend.position = "bottom", strip.text = element_text(size = 14), legend.key = element_blank())
  
  
  pdf(paste0(name, "lr_contribution_",gene_id, "_", snp_id,".pdf"), 10, 7)
  print(ggp)
  dev.off()
  
  
  return(invisible(lrcm))
  
}


for(i in 1:10){
  
  lr_contribution(dcr, gene_id = genes[i, "gene_id"], snp_id = genes[i, "snp_id"], name = paste0(out_name, "dm0seeker1_cr_", i, "_"))
  
}







### Plot examples where significant for DRIMSEQ and significant for SQTLSEEKER 


resm_sign <- resm[resm$adj_pvalue_cr < 0.05, ]

resm_sign <- resm_sign[order(resm_sign$pvalue_seeker, decreasing = FALSE), ]

genes <- resm_sign[, c("gene_id_cr", "snp_id_cr", "block_id_cr")]
colnames(genes) <- c("gene_id", "snp_id", "block_id")

genes <- genes[!duplicated(genes[, c("gene_id")]), ]



for(i in 1:10){
  # i = 1
  
  plotFit(dcr, gene_id = genes[i, "gene_id"], snp_id = genes[i, "snp_id"], out_dir = paste0(out_name, "dm0seeker0_cr_", i, "_"))
  
}



#####################################
### Plot p-values for genes where proportions are within/out the observed range 
#####################################


plot_prop_out_range_pvalues <- function(d, name){
  
  gene_ids <- names(d@counts)
  
  pi_within_range_genes <- lapply(1:length(d@counts), function(g){
    # g = 1
    
    pi_within_range_blocks <- lapply(1:nrow(d@genotypes[[g]]), function(b){
      # g = 1 ; b = 1
      
      genotypes <- d@genotypes[[g]][b, ]
      
      pi <- d@fit_full[[g]][[b]]
      
      valid_genotypes <- as.numeric(colnames(pi))
      which_valid_genotypes <- which(!is.na(pi[1, ]))
      
      pi_within_range <- rep(NA, ncol(pi))
      names(pi_within_range) <- colnames(pi)
      
      for(i in which_valid_genotypes){
        # i = 1
        
        y <- d@counts[[g]][, genotypes == valid_genotypes[i] & !is.na(genotypes)]
        
        pi_init_m <- y/matrix(colSums(y), nrow = nrow(y), ncol = ncol(y), byrow = TRUE)
        
        pi_init_minmax <- apply(pi_init_m, MARGIN = 1, function(x){  
          # c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
          range(x, na.rm = TRUE)
        })
        
        pi_tmp <- pi[, as.character(valid_genotypes[i])]
        
        pi_within_range[i] <- all(pi_tmp > pi_init_minmax[1, ] & pi_tmp < pi_init_minmax[2, ])
        
        
      }
      
      return(pi_within_range)
      
    })
    
    pi_within_range_blocks <- do.call(rbind, pi_within_range_blocks)
    pi_within_range_blocks <- data.frame(gene_id = gene_ids[g], block_id = rownames(d@genotypes[[g]]), pi_within_range_blocks, stringsAsFactors = FALSE)
    
    return(pi_within_range_blocks)
    
  })
  
  
  pi_within_range_genes <- rbind.fill(pi_within_range_genes)
  pi_within_range_genes$any_out_range <- rowSums(!pi_within_range_genes[, paste0("X", 0:2)], na.rm = TRUE) > 1
  
  res <- results(d)
  
  ### Keep only unique tests (snps)
  res <- res[!duplicated(res[, c("gene_id", "block_id")]), ]
  
  
  resm <- merge(res, pi_within_range_genes, by = c("gene_id", "block_id"), sort = FALSE)
  
  tt <- table(resm$any_out_range)
  
  ggp <- ggplot(resm, aes(x = pvalue, fill = any_out_range)) +
    theme_bw() +
    xlab("p-values") +
    ylab("Frequency") +
    geom_density(alpha = 0.6) +
    ggtitle(paste0(paste0(names(tt), ": ", tt), collapse = ", "))
  
  
  pdf(paste0(out_name, "validation_prop_out_range_pvalues_", name, ".pdf"))
  print(ggp)
  dev.off()
  
  return(invisible(resm))
  
}



plot_prop_out_range_pvalues(d = dcr, name = "cr")

plot_prop_out_range_pvalues(d = dncr, name = "ncr")






























