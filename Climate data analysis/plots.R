


## Average temperature change (beta)

par(cex.lab=1.5,cex.axis=1.3,mar=c(3.5,3,3,1),mgp=c(1.7,0.5,0))
image.plot(lon,lat, beta, zlim=c(-0.27,0.27),xlab='Longitude',ylab='Latitude')
map(add=TRUE,lwd=3,col=rgb(1, 1, 1, alpha = 0.5))
map(add=TRUE,lwd=1)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4])


## Test statistic value (t test)

par(cex.lab=1.5,cex.axis=1.3,mar=c(3.5,3,3,1),mgp=c(1.7,0.5,0))
image.plot(lon,lat, tval, zlim=c(-12,12),xlab='Longitude',ylab='Latitude')
map(add=TRUE,lwd=3,col=rgb(1, 1, 1, alpha = 0.5))
map(add=TRUE,lwd=1)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4])

cols = heat.colors(10)
colors.pval = heat.colors(141,alpha=0.4)[40:140]
colors.pval[1] = cols[1]
colors.pval[2:5] = cols[3]
colors.pval[6:10] = cols[4]

## Unadjusted p value

par(cex.lab=1.5,cex.axis=1.3,mar=c(3.5,3,3,5.5),mgp=c(1.7,0.5,0))
image(lon,lat, pval,col=heat.colors(141,alpha=0.4)[40:140],zlim=c(0,1),xlab='Longitude',ylab='Latitude')
image(lon,lat, pval<0.1, asp=1,col=c(NA,cols[4]),add=TRUE)
image(lon,lat, pval<0.05, asp=1,col=c(NA,cols[3]),add=TRUE)
image(lon,lat, pval<0.01, asp=1,col=c(NA,cols[1]),add=TRUE)
map(add=TRUE,lwd=1)
contour(lon, lat, pval, levels = c(0.05,0.01,0.1), add = T, col = 'darkblue', drawlabels = F, xlim = c(-180, 180), ylim = c(-90, 90), lty = 2,lwd=0.5)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4])

image.plot(lon,lat,pval,col=colors.pval,zlim=c(0,1),xlab='Longitude',ylab='Latitude',legend.only = TRUE)


## Adjusted p value

par(cex.lab=1.5,cex.axis=1.3,mar=c(3.5,3,3,5.5),mgp=c(1.7,0.5,0))
image(lon,lat, p.bh,col=heat.colors(141,alpha=0.4)[40:140],zlim=c(0,1),xlab='Longitude',ylab='Latitude')
image(lon,lat, p.bh<0.1, asp=1,col=c(NA,cols[4]),add=TRUE)
image(lon,lat, p.bh<0.05, asp=1,col=c(NA,cols[3]),add=TRUE)
image(lon,lat, p.bh<0.01, asp=1,col=c(NA,cols[1]),add=TRUE)
map(add=TRUE,lwd=1)
contour(lon, lat, p.bh, levels = c(0.05,0.01,0.1), add = T, col = 'darkblue', drawlabels = F, xlim = c(-180, 180), ylim = c(-90, 90), lty = 2,lwd=0.5)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4])

image.plot(lon,lat,p.bh,col=colors.pval,zlim=c(0,1),xlab='Longitude',ylab='Latitude',legend.only = TRUE)


