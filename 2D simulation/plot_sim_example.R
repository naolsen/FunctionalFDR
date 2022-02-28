
### Simulated fields example
# ## Base signal: 9 cones centered at {0.25, 0.50, 0.75}^2; diameter 0.20. Five cones pointing "upwards" and four cones pointing "downwards".

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
### Sampling of the random processes

###########
simpel.wiener.sampler <- function(felt, eval, rvals, samples, range = 1) {
  if (missing(samples)) samples <- ncol(rvals)
  
  lf <- length(felt)
  
  matrix(rvals, ncol = lf^2, nrow = samples) %*% 
    dnorm(sqrt((rep(felt, each = lf)- eval[1])^2  + (rep(felt, lf)- eval[2])^2) * range) 
}

# Simulate the fields
set.seed(936807)

B <- 2
noise <- matrix(rnorm(201*201), 201*201, 1)  ## We use 201 x 201 equidistant points on [0,1]x[0,1] for the background Wiener process.
L <- 255
vals0 <- array(, c(L, L, 2))

for (i in 1:L) {
  cat(i)
  for (j in 1:L) {
#    val <- simpel.wiener.sampler(seq(0,1, 0.005), c(j / (L+1), i / (L+1)), rvals = noise, range = 5, samples = 1) / 30
    
    vals0[i,j, ] <- 
      simpel.wiener.sampler(seq(0,1, 0.005), c(i / (L+1), j / (L+1)), rvals = noise, range = 5, samples = 2) / 30
}}



#######
# Plot
library(fields)


#pdf("..", width = 11, height = 6)

par(mfrow = c(2,4))

siglayers1 <- array(, dim = c(L,L,4))
siglayers1[,,1] <- vals0[,,1]
siglayers1[,,2] <- vals0[,,1] + 0.5*signal(1:255/ 256)
siglayers1[,,3] <- vals0[,,1] + 1.0*signal(1:255/ 256)
siglayers1[,,4] <- vals0[,,1] + 2.0*signal(1:255/ 256)


siglayers2 <- array(, dim = c(L,L,4))
siglayers2[,,1] <- vals0[,,2]
siglayers2[,,2] <- vals0[,,2] + 0.5*signal(1:255/ 256)
siglayers2[,,3] <- vals0[,,2] + 1.0*signal(1:255/ 256)
siglayers2[,,4] <- vals0[,,2] + 2.0*signal(1:255/ 256)

rang <- range(c(siglayers1, siglayers2))

image.plot(siglayers1[,,1], nlevel = 128, xaxt = "n", yaxt = "n", xlab =  0,
           zlim = rang)#, 

image.plot(siglayers1[,,2], nlevel = 128, xaxt = "n", yaxt = "n", xlab =  0.5,
           zlim = rang)#,

image.plot(siglayers1[,,3], nlevel = 128, xaxt = "n", yaxt = "n", xlab =  1,
           zlim = rang)#,

image.plot(siglayers1[,,4], nlevel = 128, xaxt = "n", yaxt = "n", xlab =  2,
           zlim = rang)#,

image.plot(siglayers2[,,1], nlevel = 128, xaxt = "n", yaxt = "n", xlab =  0,
           zlim = rang)#,

image.plot(siglayers2[,,2], nlevel = 128, xaxt = "n", yaxt = "n", xlab =  0.5,
           zlim = rang)#,

image.plot(siglayers2[,,3], nlevel = 128, xaxt = "n", yaxt = "n", xlab =  1,
           zlim = rang)#,

image.plot(siglayers2[,,4], nlevel = 128, xaxt = "n", yaxt = "n", xlab =  2,
           zlim = rang)#,

#dev.off()

