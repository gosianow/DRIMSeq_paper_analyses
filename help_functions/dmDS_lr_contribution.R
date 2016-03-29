##########################################################################
### Function for the contribution of each transcript to the LR test
### Works when d is a dmDStest object

library(reshape2)
library(ggplot2)


lik_contribution <- function(y, pi, gamma0){
  
  l <- rowSums(lgamma(y + pi * gamma0) - lgamma(pi * gamma0), na.rm = TRUE)
  
  return(l)
  
}


lr_contribution <- function(d, gene_id, name = NULL){
  
  gamma0 <- slot(d, d@dispersion)
  
  if(length(gamma0) > 1)
    gamma0 <- gamma0[gene_id]
  
  ## Null
  y <- d@counts[[gene_id]]
  pi <- d@fit_null[[gene_id]][, "null"]
  
  likn <- lik_contribution(y, pi, gamma0)
  
  ### Check if calculation is correct
  N <- ncol(y)
  S <- colSums(y)
  l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + sum(likn, na.rm = TRUE)
  ## Should be 0
  d@fit_null@metadata[gene_id, "lik"] - l
  
  
  ### Full
  groups <- levels(d@samples$group)
  
  likf <- sapply(1:length(groups), function(ii){
    # ii = 1
    
    group <- groups[ii]
    y <- d@counts[[gene_id]][, d@samples$group == group]
    pi <- d@fit_full[[gene_id]][, group]
    
    likfg <- lik_contribution(y, pi, gamma0)
    
    ### Check if calculation is correct
    N <- ncol(y)
    S <- colSums(y)
    l <- N * lgamma(gamma0) - sum(lgamma(S + gamma0)) + sum(likfg, na.rm = TRUE)
    ## Should be 0
    d@fit_full@metadata[gene_id, paste0(group)] - l 
    
    return(likfg)
    
  })
  
  ### LR contributions
  lrc <- 2 * (rowSums(likf, na.rm = TRUE) - likn)
  
  lrc <- data.frame(t(lrc))
  
  lrc$gene_id <- gene_id
  
  lrcm <- melt(lrc, id.vars = "gene_id", variable.name = "feature", value.name = "lr_contrib")
  
  lrc_tmp <- lrcm
  
  lr_tmp <- sum(lrc_tmp$lr_contrib, na.rm = TRUE)
  
  ggp <- ggplot(lrc_tmp, aes(x = feature, y = lr_contrib)) +
    geom_bar(stat = "identity") +
    ggtitle(paste0("LR = ", lr_tmp)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16), axis.text.y = element_text(size = 16), axis.title.x = element_blank(), axis.title.y = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14), legend.title = element_blank(), legend.position = "bottom", strip.text = element_text(size = 14), legend.key = element_blank())
  
  
  pdf(paste0(name, "lr_contribution_",gene_id, ".pdf"), 14, 7)
  print(ggp)
  dev.off()
  
  
  return(invisible(lrcm))
  
}


