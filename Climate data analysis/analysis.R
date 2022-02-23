
load('Data.RData')

## longitudes and latitudes
lon <- -179.5:179.5
lat <- -89.5:89.5

# Weights
weights <- matrix(, 180, 360)
for (i in 1:90) {
  weights[181-i,] <- weights[i,] <- sin((i-0.5)/180*pi)
}
weights <- t(weights/sum(weights))


## data and linear regression
d2 <- matrix(data[,16], nrow = 25)
L <- lm(d2 ~ c(1:25))


## yearly change
beta <- matrix(coef(L)[2,], 360, 180)

## t-value
tval <- matrix(coef(L)[2,] / sqrt(colSums(L$residuals^2)/23/1300),
               360, 180)

## unadjusted p-value
pval <- pt(tval, df = 23, lower.tail = FALSE)


## Adjustment procedure
p.order <- order(pval)
p.bh <- matrix(1:64800, 360, 180)

wcum <- cumsum(t(weights[p.order]))
for (i in 1:64800) {
  p.bh[[i]] <- pval[p.order[i]] / wcum[i]
}
for (i in 64799:1) {
  p.bh[i] <- min(p.bh[i:(i+1)])
}
p.bh[p.order] <- p.bh 
## p.bh: adjusted p-values


#### Results
## Percentage of significance, unadjusted p-values
sum((pval < 0.10)* weights)
sum((pval < 0.05)* weights)
sum((pval < 0.01)* weights) 
sum((pval < 0.001)* weights)


## Percentage of significance, adjusted p-values
sum((p.bh < 0.10)* weights)
sum((p.bh < 0.05)* weights)
sum((p.bh < 0.01)* weights) 
sum((p.bh < 0.001)* weights)




