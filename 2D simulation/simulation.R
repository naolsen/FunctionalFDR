

#### Simulation on FDR in functional/spatial setting
# Niels Olsen, juni 2017


# --- Note: 
# This file is intended for execution in parallel (using doParallel) and takes a long time to run.
# Can possibly be implemented in a muuch faster way.


## Set up backend
library(foreach)
library(doParallel)
registerDoParallel(cores = 24)



## Base signal: 9 cones centered at {0.25, 0.50, 0.75}^2; diameter 0.20. Five cones pointing "upwards" and four cones pointing "downwards".

# Evaluates the base signal at specified time points
signal <- function(time) {
  l <- length(time)
  ud <- matrix(0, l,l)
  
  w1 <- which(time > 0.15 &  time < 0.35)
  w2 <- which(time > 0.40 &  time < 0.60)
  w3 <- which(time > 0.65 &  time < 0.85)
  
  ud[w1,w1] <- 
    10 * pmax.int(0.1 - sqrt(outer( (0.25 - time[w1])^2, (0.25 - time[w1])^2, "+")), 0)
  
  ud[w2,w1] <- t(ud[w1,w2] <- 
                   - 10 * pmax.int(0.1 - sqrt(outer( (0.25 - time[w1])^2, (0.50 - time[w2])^2, "+")), 0))
  
  ud[w3,w1] <- t(ud[w1,w3] <-  
                   10 * pmax.int(0.1 - sqrt(outer( (0.25 - time[w1])^2, (0.75 - time[w3])^2, "+")), 0))
  
  ud[w2,w2] <- 
    10 * pmax.int(0.1 - sqrt(outer( (0.50 - time[w2])^2, (0.50 - time[w2])^2, "+")), 0)
  
  ud[w3,w2] <- t(ud[w2,w3] <- 
                   - 10 * pmax.int(0.1 - sqrt(outer( (0.50 - time[w2])^2, (0.75 - time[w3])^2, "+")), 0))
  
  ud[w3,w3] <- 
    10 * pmax.int(0.1 - sqrt(outer( (0.75 - time[w3])^2, (0.75 - time[w3])^2, "+")), 0)
  
  ud
}

contour(signal(1:500 / 500))



### Sampling of the random processes


###########
simpel.wiener.sampler <- function(felt, eval, rvals, samples, range = 1) {
  if (missing(samples)) samples <- ncol(rvals)
  
  lf <- length(felt)
  
  matrix(rvals, ncol = lf^2, nrow = samples) %*% 
    dnorm(sqrt((rep(felt, each = lf)- eval[1])^2  + (rep(felt, lf)- eval[2])^2) * range) 
}




######
L <- 255 # Grid size
M <- 1000 # No of replicates

scenarios <- list(c(2, 20), c(2, 10), c(2, 40), c(1, 20), c(0.5, 20))


results <- list()

for (m in 1:length(scenarios)) {
  
  signal.size <- scenarios[[m]][[1]] ## Size of signal
  fixed.signal <- signal(1:255/ 256)
  sig0 <- which(signal(1:255/ 256) != 0) # Where the null hypothesis gets rejected
  
  B <-  scenarios[[m]][2] ## No. of samples
  
  plist <- 
    foreach(r = 1:M, .combine = cbind) %dopar% {
      noise <- matrix(rnorm(201*201*B), 201*201, B)  ## We use 201 x 201 equidistant points on [0,1]x[0,1] for the background Wiener process.
      
      pvals <- matrix(, L, L)
      
      for (i in 1:L) {
        cat(i)
        for (j in 1:L) {
          val <- simpel.wiener.sampler(seq(0,1, 0.005), c(j / (L+1), i / (L+1)), rvals = noise, range = 5, samples = B) / 30 + signal.size*fixed.signal[i,j]
          pvals[i,j] <- t.test(val)$p.val
        }}
      pvals
    }
  dim(plist) <- c(L, L, M)
  
  ## Scenario 
  
  a.levels <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.10)
  
  
  FDR.u <- FDR <-  power <- signif <- matrix(NA, M,7)
  
  pvals <- plist
  
  for (i in 1:M) {
    pval.mat <- pvals[,,i]
    
    pval.stor <- p.adjust(pval.mat, method = "BH")
    
    for (l in 1:7) {
      
      ## Power / sensitivity
      power[i,l] <- sum((pval.stor <= a.levels[l])[sig0]) 
      
      ## Significance / type I errors
      signif[i,l] <-  sum((pval.stor <= a.levels[l])[-sig0]) 
      
      ## FDR
      FDR[i,l] <- sum((pval.stor <= a.levels[l])[-sig0]) / sum((pval.stor <= a.levels[l]))
      
      ## unadjusted FDR
      FDR.u[i,l] <- sum((pval.mat <= a.levels[l])[-sig0]) / sum((pval.mat <= a.levels[l]))
    }}
  FDR[is.nan(FDR)] <- 0
  
  results[[m]] <- list(fdr = FDR, power = power, fdr.u = FDR.u, signific = signif)
}


save(results, file = "FDR_spike_sim_pvals_list.RData")



