
calculateCAT <- function(results1, results2, by = 1, maxrank = min(nrow(results1), nrow(results1)), FDR = 0.05){
  
  
  results1 <- results1[order(results1$pvalue, decreasing = FALSE), , drop = FALSE]
  results2 <- results2[order(results2$pvalue, decreasing = FALSE), , drop = FALSE]
  
  results1 <- results1[complete.cases(results1), , drop = FALSE]
  results2 <- results2[complete.cases(results2), , drop = FALSE]
  
  x1 <- min(which(!results1$adj_pvalue < FDR))
  x2 <- min(which(!results2$adj_pvalue < FDR))
  
  
  top_ds_genes <- sort(unique(c(seq(1, maxrank, by = by), x1, x2)), decreasing = FALSE)
  overlap <- rep(0, length(top_ds_genes))
  
  overpal_nr_so_far <- 0
  overpal_elem_so_far <- NULL
  
  for(i in 1:length(top_ds_genes)){
    # i = 1000
    
    set1 <- results1[1:top_ds_genes[i], "gene_id"]
    set2 <- results2[1:top_ds_genes[i], "gene_id"]
    
    overlap_elem <- intersect(set1, set2)
    
    overlap[i] <- length(overlap_elem) / top_ds_genes[i]
    
  }
  
  overlaps <- data.frame(top_ds_genes = top_ds_genes, overlap = overlap)
  
  mm <- match(c(x1, x2), top_ds_genes)
  
  X <- data.frame(top_ds_genes = c(x1, x2), overlap = overlaps$overlap[mm], results = c("results1", "results2"), stringsAsFactors = FALSE)
  
  return(list(overlaps = overlaps, X = X))
  
}




plotCAT <- function(data_CAT, metadata, plot_var, facet_var = NULL, plot_colors = NULL, plotx = TRUE, reference_color = NULL){
  
  stopifnot(length(data_CAT) == nrow(metadata))
  
  stopifnot(all(sapply(data_CAT, function(x){
    
    is.list(x) && length(x) == 2 && all(names(x) == c("overlaps", "X"))
    
    })))
  
  stopifnot(plot_var %in% colnames(metadata))
  if(!is.null(facet_var))
    stopifnot(facet_var %in% colnames(metadata))
  
  metadata <- metadata[, c(plot_var, facet_var), drop = FALSE]
  metadata[, plot_var] <- factor(metadata[, plot_var])
  
  if(!is.null(plot_colors))
    stopifnot(nlevels(metadata[, plot_var]) == length(plot_colors))
  
  
  overlapslist <- lapply(1:length(data_CAT), function(r){
    # r = 1
    
    overlaps <- data_CAT[[r]]$overlaps
    
    overlaps <- cbind(overlaps, metadata[rep(r, nrow(overlaps)), , drop = FALSE])
    
  })
  
  Xlist <- lapply(1:length(data_CAT), function(r){
    # r = 1
    
    X <- data_CAT[[r]]$X
    
    X <- cbind(X, metadata[rep(r, nrow(X)), , drop = FALSE])
    
  })
  
  overlaps <- rbind.fill(overlapslist)
  X <- rbind.fill(Xlist)
  
  
  overlaps$plot_var <- overlaps[, plot_var]
  
  X1 <- subset(X, results == "results1", select = colnames(X) != "results") # reference
  X2 <- subset(X, results == "results2", select = colnames(X) != "results")
  
  X2$plot_var <- X2[, plot_var]
  
  ggp <- ggplot(data = overlaps, aes(x = top_ds_genes, y = overlap, group = plot_var, colour = plot_var)) +
    theme_bw() +
    geom_line(size = 2, na.rm=TRUE) +
    xlab("Number of top ranked genes") +
    ylab("Percentage overlap") +
    theme(axis.text=element_text(size = 16), axis.title = element_text(size = 18, face = "bold"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 12), strip.text = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size = 1.5, shape = NA), ncol = 3)) 
  
  if(plotx)
    ggp <- ggp + geom_point(data = X2, aes(x = top_ds_genes, y = overlap, colour = plot_var), size = 10, shape = "X", na.rm=TRUE) 
  
  if(plotx)
    ggp <- ggp + geom_point(data = X1, aes(x = top_ds_genes, y = overlap), colour = reference_color, size = 8, shape = "X", na.rm=TRUE) 
  
  if(!is.null(plot_colors))
    ggp <- ggp + scale_color_manual(values = plot_colors)
  
  if(length(facet_var) == 1)
    ggp <- ggp + facet_wrap(reformulate(facet_var[1]))
  
  if(length(facet_var) == 2)
    ggp <- ggp + facet_grid(reformulate(facet_var[1], facet_var[2]))
  
  
  return(ggp)
  
  
}



