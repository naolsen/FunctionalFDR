

## Base signal: 9 cones centered at {0.25, 0.50, 0.75}^2; diameter 0.20. Five cones pointing "upwards" and four cones pointing "downwards".

### Evaluates the base signal at specified time points

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



## Analysis

load("FDR_spike_sim_pvals_list.RData")


fixed.signal <- signal(1:255/ 256)
sig0 <- which(signal(1:255/ 256) != 0)
  
## Threshold levels
a.levels <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.10)



#### -------
  
  ## Scenario 1
  
  FDR <-  power <- signif <- matrix(NA, 250,7)
  FDR.u <- matrix(NA, 250,7)


  FDR <- results[[1]]$fdr
  FDR.u <- results[[1]]$fdr.u
  signif <- results[[1]]$signific
  power <- results[[1]]$power

  FDR[is.nan(FDR)] <- 0
  FDR.u[is.nan(FDR.u)] <- 0

  colMeans(power / length(sig0))
  colMeans(signif / (255*255 - length(sig0)))
  colMeans(FDR)
  colMeans(FDR.u)
  
  eksp1 <- cbind(power = colMeans(power / length(sig0)), 
                 significance = colMeans(signif / (255*255 - length(sig0))), 
                 FDR = colMeans(FDR),FDR.u = colMeans(FDR.u), deparse.level = 1)
  

  ## Scenario 2 -----
  
  FDR <- results[[2]]$fdr
  FDR.u <- results[[2]]$fdr.u
  signif <- results[[2]]$signific
  power <- results[[2]]$power
  
  FDR[is.nan(FDR)] <- 0
  FDR.u[is.nan(FDR.u)] <- 0
  
  colMeans(power / length(sig0))
  colMeans(signif / (255*255 - length(sig0)))
  colMeans(FDR) 
  colMeans(FDR.u)
  
  
  eksp2 <- cbind(power = colMeans(power / length(sig0)), 
                 significance = colMeans(signif / (255*255 - length(sig0))), 
                 FDR = colMeans(FDR),FDR.u = colMeans(FDR.u), deparse.level = 1)
  
  
  ## Scenario 3 -----
  
  FDR <- results[[3]]$fdr
  FDR.u <- results[[3]]$fdr.u
  signif <- results[[3]]$signific
  power <- results[[3]]$power
  
  FDR[is.nan(FDR)] <- 0
  FDR.u[is.nan(FDR.u)] <- 0
  
  colMeans(power / length(sig0))
  colMeans(signif / (255*255 - length(sig0)))
  colMeans(FDR) 
  colMeans(FDR.u)
  
  
  eksp3 <- cbind(power = colMeans(power / length(sig0)), 
                 significance = colMeans(signif / (255*255 - length(sig0))), 
                 FDR = colMeans(FDR),FDR.u = colMeans(FDR.u), deparse.level = 1)
  
  ## Scenario 4 -----
  
  FDR <- results[[4]]$fdr
  FDR.u <- results[[4]]$fdr.u
  signif <- results[[4]]$signific
  power <- results[[4]]$power
  
  FDR[is.nan(FDR)] <- 0
  FDR.u[is.nan(FDR.u)] <- 0
  
  colMeans(power / length(sig0))
  colMeans(signif / (255*255 - length(sig0)))
  colMeans(FDR) 
  colMeans(FDR.u)
  
  
  eksp4 <- cbind(power = colMeans(power / length(sig0)), 
                 significance = colMeans(signif / (255*255 - length(sig0))), 
                 FDR = colMeans(FDR),FDR.u = colMeans(FDR.u), deparse.level = 1)
  
  ## Scenario 5 -----
  FDR <- results[[5]]$fdr
  FDR.u <- results[[5]]$fdr.u
  signif <- results[[5]]$signific
  power <- results[[5]]$power
  
  FDR[is.nan(FDR)] <- 0
  FDR.u[is.nan(FDR.u)] <- 0
  
  colMeans(power / length(sig0))
  colMeans(signif / (255*255 - length(sig0)))
  colMeans(FDR) 
  colMeans(FDR.u)
  
  
  eksp5 <- cbind(power = colMeans(power / length(sig0)), 
                 significance = colMeans(signif / (255*255 - length(sig0))), 
                 FDR = colMeans(FDR),FDR.u = colMeans(FDR.u), deparse.level = 1)

#save(eksp1, eksp2, eksp3, eksp4, eksp5, file ="sim_results_1000.RData")

  
#### Bootstraping confidence intervals.
# Scenario 3 example -----
  
FDR <- results[[3]]$fdr

vals <- numeric(10000)  
for (j in 1:10000)  
vals[j] <-  mean(FDR[sample(1000, replace = T), 6])  
hist(vals)

sort(vals)[c(250, 9750)]

## Bootstrap values, 10000 resamples
out <- array(, dim = c(5, 7))
for (i in 1:5) for (j in 1:7) {
  
  FDR <- results[[i]]$fdr
  
  vals <- numeric(10000)  
  for (k in 1:10000)  
    vals[k] <-  mean(FDR[sample(1000, replace = T), j])  
  
  # CI or std. error
  
  #out[,i,j] <- sort(vals)[c(250, 9750)] 
  out[i,j] <- sd(vals)
}
out


  
  
  
  
  
  
  
