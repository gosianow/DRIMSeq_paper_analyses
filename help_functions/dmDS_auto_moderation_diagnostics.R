library(ggplot2)
library(DRIMSeq)

#' Function to produce plots about automatic moderation 
#' 
#' @param x dmDSdispersion object

dmDS_auto_moderation_diagnostics <- function(x, common_disp, out_dir_tmp, BPPARAM){
  

  mean_expression = TRUE
  common_dispersion = FALSE
  genewise_dispersion = TRUE
  disp_adjust = TRUE
  disp_mode = "grid"
  disp_interval = c(0, 1e+05)
  disp_tol = 1e-08
  disp_init = common_disp
  disp_init_weirMoM = TRUE
  disp_grid_length = 21
  disp_grid_range = c(-10, 10)
  disp_moderation = "none"
  disp_prior_df = 0
  disp_span = 0.1
  prop_mode = "constrOptimG"
  prop_tol = 1e-12
  verbose = TRUE
  
  
  ###### Inside the dmDispersion function
  
  mean_expression <- DRIMSeq:::dm_estimateMeanExpression(counts = x@counts, verbose = verbose, BPPARAM = SerialParam())
  
  
  counts = x@counts
  samples = x@samples
  
  
  
  ###### Inside the dmDS_estimateTagwiseDispersion function
  
  
  inds <- 1:length(counts)
  
  group <- samples$group
  ngroups <- nlevels(group)
  lgroups <- levels(group)
  nlibs <- length(group)
  
  igroups <- lapply(lgroups, function(gr){which(group == gr)})
  names(igroups) <- lgroups
  
  
  
  ###### Switch to grid
  
  splinePts <- seq(from = disp_grid_range[1], to = disp_grid_range[2], length = disp_grid_length)
  
  disp_grid_length <- length(splinePts)
  
  splineDisp <- disp_init * 2^splinePts
  
  
  ### calculate the likelihood for each gene at the spline dispersion points
  seq_disp_grid_length <- seq(disp_grid_length)
  
  message("Doing loglik")
  
  loglikL <- BiocParallel::bplapply(inds, function(g){
    # g = 1237
    # print(g)
    
    ll <- numeric(disp_grid_length)
    
    for(i in seq_disp_grid_length){
      # i = 1
      
      out <- DRIMSeq:::dm_profileLikTagwise(gamma0 = splineDisp[i], y = counts[[g]], ngroups = ngroups, lgroups = lgroups, igroups = igroups, disp_adjust = disp_adjust, prop_mode = prop_mode, prop_tol = prop_tol, verbose = verbose)
      
      if(is.na(out)){
        ll <- rep(NA, disp_grid_length)
        break
      }
      
      ll[i] <- out
      
    }
    
    return(ll)
    
  }, BPPARAM = BPPARAM)
  
  message("Loglik done!")
  
  loglik <- do.call(rbind, loglikL)
  
  not_nas <- complete.cases(loglik)        
  
  loglik <- loglik[not_nas, , drop = FALSE]
  mean_expression <- mean_expression[not_nas]
  
  
  loglik_orig <- loglik
  
  
  
  
  ### Check where the grid is maximized 
  grid_max <- apply(loglik, 1, which.max)
  lik_max <- apply(loglik, 1, max)
  
  ### Calculate the span of loglik
  loglik_span <- apply(loglik, 1, function(x){max(x) - min(x)})
  
  
  ### In the calculation of moderation, do not take into account genes that have dispersion on the top and bottom boundry of the grid (skipp 4 last grid points and 1 first grid point)
  not_boundry <- grid_max < (disp_grid_length - 3) & grid_max > 1
  boundry_last <- grid_max == disp_grid_length
  table(boundry_last)
  
  
  
  ############################################################################
  
  
  prefix_mod <- paste0(out_dir_tmp, "none_")
  
  
  # Plot concentration versus mean with grid points

  ggp <- plotDispersion(x)
  
  ggp <- ggp + 
    geom_abline(intercept = log10(splineDisp), slope = 0, linetype = 3, color = "orange")
  
  pdf(paste0(prefix_mod, "moderation_none_dispersion_vs_mean_grid.pdf"))
  print(ggp)
  dev.off()
  
  png(paste0(prefix_mod, "moderation_none_dispersion_vs_mean_grid.png"))
  print(ggp)
  dev.off()
  
  
  # Plot lik span versus mean expression only for the boundry genes
  df_loglik_span_mean <- data.frame(loglik_span = loglik_span[boundry_last], mean_expression = mean_expression[boundry_last])
  
  ggp <- ggplot(df_loglik_span_mean, aes(x = log10(mean_expression), y = loglik_span)) +
    geom_point(alpha = 1, size = 1) +
    coord_cartesian(xlim = log10(range(mean_expression))) 
  
  pdf(paste0(prefix_mod, "liks_span_versus_mean_expression_boundry.pdf"))
  print(ggp)
  dev.off()
  
  
  # Plot loglikelihood for one gene that is OK and one that is on the boundry of the grid
  
  
  loglik_tmp_boundry <- loglik[which(boundry_last)[1], ]
  loglik_tmp_notboundry <- loglik[which(not_boundry)[1], ]
  

  df_liks <- data.frame(splinePts = splinePts, loglik = loglik_tmp_boundry)

  ggp <- ggplot(df_liks, aes(x = splinePts, y = -log10(-loglik))) +
    geom_line(size = 3) +
    ggtitle(paste0("Loglik span ", max(df_liks$loglik) - min(df_liks$loglik))) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 14), legend.position = "bottom") 
  
  pdf(paste0(prefix_mod, "liks_boundry.pdf"))
  print(ggp)
  dev.off()
  
  
  
  df_liks <- data.frame(splinePts = splinePts, loglik = loglik_tmp_notboundry)
  
  ggp <- ggplot(df_liks, aes(x = splinePts, y = -log10(-loglik))) +
    geom_line(size = 3) +
    ggtitle(paste0("Loglik span ", max(df_liks$loglik) - min(df_liks$loglik))) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 14), legend.position = "bottom") 
  
  pdf(paste0(prefix_mod, "liks_notboundry.pdf"))
  print(ggp)
  dev.off()
  
  
  
  
  ############################################################################
  ### Auto moderation for adjustment to common dispersion
  
  prefix_mod <- paste0(out_dir_tmp, "common_")
  
  
  
  ## Calculate common moderation
  if(sum(not_boundry) == length(not_boundry)){
    moderation <- colMeans(loglik)
  }else{
    moderation <- colMeans(loglik[not_boundry, , drop = FALSE])
  }
  
  common_moderation_span <- moderation_span <- max(moderation) - min(moderation)
  
  
  
  
  # Plot boxplot of likelihood span for boundry and not boundry genes
  
  df_loglik_span <- data.frame(loglik_span = loglik_span, boundry = boundry_last)
  
  ggp <- ggplot(df_loglik_span[df_loglik_span$loglik_span > 1e-5, ], aes(x = boundry, y = log10(loglik_span), fill = boundry)) +
    geom_boxplot(outlier.size = 0) +
    geom_point(position = position_jitter(width = 0.5), alpha = 0.2) +
    geom_abline(intercept = log10(moderation_span), slope = 0, linetype = 2) 
  
  pdf(paste0(prefix_mod, "liks_span_log.pdf"))
  print(ggp)
  dev.off()
  
  
  ggp <- ggplot(df_loglik_span, aes(x = boundry, y = loglik_span, fill = boundry)) +
    geom_boxplot(outlier.size = 0) +
    geom_point(position = position_jitter(width = 0.5), alpha = 0.2) +
    geom_abline(intercept = moderation_span, slope = 0, linetype = 2) +
    coord_cartesian(ylim = c(0, 2*moderation_span))
  
  pdf(paste0(prefix_mod, "liks_span.pdf"))
  print(ggp)
  dev.off()
  
  
  
  
  
  # Plot lik span versus mean expression 
  
  df_loglik_span_mean <- data.frame(loglik_span = loglik_span, mean_expression = mean_expression, boundry = boundry_last)
  
  ggp <- ggplot(df_loglik_span_mean[df_loglik_span_mean$loglik_span > 1e-5, ], aes(x = log10(mean_expression), y = log10(loglik_span), color = boundry)) +
    geom_point(alpha = 1, size = 0.5) +
    geom_abline(intercept = log10(moderation_span), slope = 0, linetype = 2)
  
  pdf(paste0(prefix_mod, "liks_span_versus_mean_expression_log.pdf"))
  print(ggp)
  dev.off()
  
  
  ggp <- ggplot(df_loglik_span_mean, aes(x = log10(mean_expression), y = loglik_span, color = boundry)) +
    geom_point(alpha = 1, size = 0.5) +
    geom_abline(intercept = moderation_span, slope = 0, linetype = 2) +
    coord_cartesian(ylim = c(0, 2*moderation_span))
  
  pdf(paste0(prefix_mod, "liks_span_versus_mean_expression.pdf"))
  print(ggp)
  dev.off()
  
  
  
  
  
  
  
  ### Calculate the ratio between moderation lik span and lik span of boundry genes
  
  loglik_span_boundry <- df_loglik_span_mean[boundry_last, "loglik_span"]
  
  span_ratio <- moderation_span / loglik_span_boundry
  
  priorN <- quantile(1/span_ratio, 0.5)
  
  message(paste0("! Using ", round(priorN, 2), " as a shrinkage factor !"))
  
  write.table(priorN, paste0(prefix_mod, "priorn_median.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  
  
  ### Plot priorN versus mean expression
  
  df_priorN <- data.frame(priorN = 1/span_ratio, mean_expression = df_loglik_span_mean[boundry_last, "mean_expression"])
  
  
  ggp <- ggplot(df_priorN, aes(x = log10(mean_expression), y = priorN)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = priorN, slope = 0, linetype = 2) +
    coord_cartesian(xlim = log10(range(df_loglik_span_mean$mean_expression))) +
    ggtitle(paste0("priorN = ", round(priorN, 4))) 
  
  pdf(paste0(prefix_mod, "priorn_versus_mean_expression.pdf"))
  print(ggp)
  dev.off()
  
  
  ggp <- ggplot(df_priorN, aes(x = log10(mean_expression), y = log10(priorN))) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = log10(priorN), slope = 0, linetype = 2) +
    coord_cartesian(xlim = log10(range(df_loglik_span_mean$mean_expression))) +
    ggtitle(paste0("priorN = ", round(priorN, 4))) 
  
  pdf(paste0(prefix_mod, "priorn_versus_mean_expression_log.pdf"))
  print(ggp)
  dev.off()
  
  
  
  ### Fit loess to log of priorN and log of mean expression - nice results
  
  df_priorN_loglog <- data.frame(priorN = log10(1/span_ratio), mean_expression = log10(df_loglik_span_mean[boundry_last, "mean_expression"]))
  
  priorN_loess_loglog <- loess(priorN ~ mean_expression, df_priorN_loglog, control = loess.control(surface = "direct"))
  priorN_predict_loglog <- predict(priorN_loess_loglog, data.frame(mean_expression = log10(mean_expression)), se = FALSE)
  
  df_priorN_gene_loglog <- data.frame(mean_expression = mean_expression, priorN = 10 ^ priorN_predict_loglog)
  
  
  ggp <- ggplot(df_priorN, aes(x = log10(mean_expression), y = priorN)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = priorN, slope = 0, linetype = 2) +
    coord_cartesian(xlim = log10(range(df_loglik_span_mean$mean_expression))) +
    ggtitle(paste0("priorN = ", round(priorN, 4))) +
    geom_point(data = df_priorN_gene_loglog, aes(x = log10(mean_expression), y = priorN), color = "blue", size = 1)+
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 14), legend.position = "bottom") 
  
  
  pdf(paste0(prefix_mod, "priorn_versus_mean_expression_loess.pdf"))
  print(ggp)
  dev.off()
  
  png(paste0(prefix_mod, "priorn_versus_mean_expression_loess.png"))
  print(ggp)
  dev.off()
  
  
  
  
  
  
  
  ############################################################################
  ### Auto moderation for adjustment to trended dispersion
  
  loglik <- loglik_orig
  
  prefix_mod <- paste0(out_dir_tmp, "trended_")
  
  
  
  moderation <- DRIMSeq:::dm_movingAverageByCol(loglik = loglik, mean_expression = mean_expression, not_boundry = not_boundry, disp_span = disp_span)
  
  moderation_span <- apply(moderation, 1, function(x){max(x) - min(x)})
  
  
  
  # Plot lik span versus mean expression 
  
  df_loglik_span_mean <- data.frame(loglik_span = loglik_span, mean_expression = mean_expression, boundry = boundry_last)
  
  df_moderation_span <- data.frame(moderation_span = moderation_span, mean_expression = mean_expression, boundry = boundry_last)
  
  ggp <- ggplot(df_loglik_span_mean[df_loglik_span_mean$loglik_span > 1e-5, ], aes(x = log10(mean_expression), y = log10(loglik_span), color = boundry)) +
    geom_point(alpha = 1, size = 0.5) +
    geom_point(data = df_moderation_span, aes(x = log10(mean_expression), y = log10(moderation_span)), color = "black", size = 0.5) +
    geom_abline(intercept = log10(common_moderation_span), slope = 0, linetype = 3, size = 2) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 14), legend.position = "bottom") 
  

  pdf(paste0(prefix_mod, "liks_span_versus_mean_expression_log.pdf"))
  print(ggp)
  dev.off()
  
  png(paste0(prefix_mod, "liks_span_versus_mean_expression_log.png"))
  print(ggp)
  dev.off()
  
  
  ggp <- ggplot(df_loglik_span_mean, aes(x = log10(mean_expression), y = loglik_span, color = boundry)) +
    geom_point(alpha = 1, size = 0.5) +
    geom_point(data = df_moderation_span, aes(x = log10(mean_expression), y = moderation_span), color = "black", size = 0.5) +
    coord_cartesian(ylim = c(0, 2*max(moderation_span)))
  
  pdf(paste0(prefix_mod, "liks_span_versus_mean_expression.pdf"))
  print(ggp)
  dev.off()
  
  
  
  ### Calculate the ratio between moderation lik span and lik span of boundry genes
  
  loglik_span_boundry <- df_loglik_span_mean[boundry_last, "loglik_span"]
  moderation_span_boundry <- df_moderation_span[boundry_last, "moderation_span"]
  
  
  span_ratio <- moderation_span_boundry / loglik_span_boundry
  
  priorN <- quantile(1/span_ratio, 0.5)
  
  message(paste0("! Using ", round(priorN, 2), " as a shrinkage factor !"))
  
  write.table(priorN, paste0(prefix_mod, "priorn_median.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  
  
  
  ### Plot priorN versus mean expression
  
  df_priorN <- data.frame(priorN = 1/span_ratio, mean_expression = df_loglik_span_mean[boundry_last, "mean_expression"])
  
  
  ggp <- ggplot(df_priorN, aes(x = log10(mean_expression), y = priorN)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = priorN, slope = 0, linetype = 2) +
    coord_cartesian(xlim = log10(range(df_loglik_span_mean$mean_expression))) +
    ggtitle(paste0("priorN = ", round(priorN, 4))) 
  
  pdf(paste0(prefix_mod, "priorn_versus_mean_expression.pdf"))
  print(ggp)
  dev.off()
  
  
  ggp <- ggplot(df_priorN, aes(x = log10(mean_expression), y = log10(priorN))) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = log10(priorN), slope = 0, linetype = 2) +
    coord_cartesian(xlim = log10(range(df_loglik_span_mean$mean_expression))) +
    ggtitle(paste0("priorN = ", round(priorN, 4)))
  
  
  pdf(paste0(prefix_mod, "priorn_versus_mean_expression_log.pdf"))
  print(ggp)
  dev.off()
  
  
  
  ### Fit loess to log of priorN and log of mean expression - nice results
  
  df_priorN_loglog <- data.frame(priorN = log10(1/span_ratio), mean_expression = log10(df_loglik_span_mean[boundry_last, "mean_expression"]))
  
  priorN_loess_loglog <- loess(priorN ~ mean_expression, df_priorN_loglog, control = loess.control(surface = "direct"))
  priorN_predict_loglog <- predict(priorN_loess_loglog, data.frame(mean_expression = log10(mean_expression)), se = FALSE)
  
  df_priorN_gene_loglog <- data.frame(mean_expression = mean_expression, priorN = 10 ^ priorN_predict_loglog)
  
  
  
  ggp <- ggplot(df_priorN, aes(x = log10(mean_expression), y = priorN)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = priorN, slope = 0, linetype = 2) +
    coord_cartesian(xlim = log10(range(df_loglik_span_mean$mean_expression))) +
    ggtitle(paste0("priorN = ", round(priorN, 4))) +
    geom_point(data = df_priorN_gene_loglog, aes(x = log10(mean_expression), y = priorN), color = "blue", size = 1) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(size = 14), legend.position = "bottom") 
  
  pdf(paste0(prefix_mod, "priorn_versus_mean_expression_loess.pdf"))
  print(ggp)
  dev.off()
  
  png(paste0(prefix_mod, "priorn_versus_mean_expression_loess.png"))
  print(ggp)
  dev.off()
  
  return("All done!!! :) ")
  
  
  }










































