############################################################################
# function to simulate data from DM
# Created 28 Oct 2015
############################################################################

library(parallel)

# m = 10; n = 5; pi = c(1/3, 1/3, 1/3); g0 = 100; nm = 100; tot = "uni"; nd = 3; BPPARAM = BiocParallel::MulticoreParam(workers = 4)

dm_rdirichlet <- function(n = 1, alpha){
  x <- matrix(0,n,length(alpha))
  for(i in 1:length(alpha)) 
    x[,i] <- rgamma(n, shape=alpha[i])
  x/rowSums(x)
}


dm_simulate <- function(m = 10, n = 5, pi = c(1/3, 1/3, 1/3), g0 = 100, nm = 100, tot = "uni", nd = 0.25, mc.cores = 1){
  # m - number of genes
  # n - total number of samples/replicates
  # pi - proportions of features 
  # g0 - dispersion gamma0 (can be a vector of length m)
  # nm - total number of counts per gene (can be a vector of length m)
  # tot - method to generate total number of counts per gene: nbinom, norm, uni
  # nd - dispersion of nm if simulated from nbinom or norm (can be a vector of length m)
  
  if(length(g0 == 1))
    g0 <- rep(g0, m)
  if(length(nm == 1))
    nm <- rep(nm, m)
  if(length(nd == 1))
    nd <- rep(nd, m)
  
  sim <- mclapply(1:m, function(i){
    # i=1
    
    g_dir_org <- pi * g0[i] # gamma
    
    # simulate dirichlets
    g_dir <- dm_rdirichlet(n, g_dir_org)
    
    # simulate total counts
    if(tot == "nbinom")
      t <- rnbinom(n, mu = nm[i], size = 1/nd[i]) # Variance = mu + mu^2/size
#     if(tot == "norm")
#       t <- round(rnorm(n, mean = nm[i], sd = nd[i]))  
    if(tot == "uni")
      t <- rep(nm[i], n)
    
    # simulate multinomial
    dm <- sapply(1:n, function(j) rmultinom(1, prob = g_dir[j, ], size = t[j]))    
   
    rownames(dm) <- paste0("g", rep(i, length(g_dir_org)), ":f", 1:length(g_dir_org))
    
    return(dm)
    
  }, mc.cores = mc.cores)
  
sim <- do.call(rbind, sim)

sim

}









