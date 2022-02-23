
#### This file plots the results of Simulation 1 in one large figure

layout(cbind(rbind(1:5,6:10,11:15,c(0,16,16,16,0))),heights=c(0.9,0.9,0.9,0.3))
par(mar=c(2,2.5,2,0.1),font.lab=2,cex.lab=1.2,mgp=c(1.2,0.45,0),cex.main=1.4)

## Scenario 1
load('sim_scen1_M1000_B5000.RData')
matplot(ascissa,t(y_i),ylab='',type='l',col=gray((1-x_i)*0.6+0.2),lwd=1,lty=1,main='Data',ylim=c(-2,6.5),xlab='')
title(ylab = "Scenario 1", mgp = c(1.5, 0.4, 0))

plot(d,FWER.unadj,ylim=c(0,1),main='FWER',type='l',lwd=2,lty=1,pch=1,col=1,
     ylab='',xlab='d')
lines(d,FWER.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,FWER.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
grid()
abline(h=0.05,lty=1)

plot(d,FDR.unadj,ylim=c(0,1),main='FDR',type='l',lwd=2,lty=1,pch=1,col=1,xlab='d',ylab='')
lines(d,FDR.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,FDR.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
grid()
abline(h=0.05,lty=1)

plot(d,spec.unadj,ylim=c(0,1),main='Specificity',type='l',lwd=2,lty=1,pch=1,col=1,xlab='d',ylab='')
lines(d,spec.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,spec.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
grid()

plot(d,sens.unadj,ylim=c(0,1),main='Sensitivity',type='l',lwd=2,lty=1,pch=1,col=1,xlab='d',ylab='')
lines(d,sens.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,sens.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
grid()

## Scenario 2
load('sim_scen2_M1000_B5000.RData')
matplot(ascissa,t(y_i),type='l',ylab='',col=gray((1-x_i)*0.6+0.2),lwd=1,lty=1,main='',ylim=c(-2,6.5),xlab='')
title(ylab = "Scenario 2", mgp = c(1.5, 0.4, 0))

plot(d,FWER.unadj,ylim=c(0,1),main='',type='l',lwd=2,lty=1,pch=1,col=1,
     ylab='',xlab='d')
lines(d,FWER.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,FWER.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
grid()
abline(h=0.05,lty=1)

plot(d,FDR.unadj,ylim=c(0,1),main='',type='l',lwd=2,lty=1,pch=1,col=1,xlab='d',ylab='')
lines(d,FDR.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,FDR.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
grid()
abline(h=0.05,lty=1)

plot(d,spec.unadj,ylim=c(0,1),main='',type='l',lwd=2,lty=1,pch=1,col=1,xlab='d',ylab='')
lines(d,spec.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,spec.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
abline(h=0.05,lty=1)
grid()

plot(d,sens.unadj,ylim=c(0,1),main='',type='l',lwd=2,lty=1,pch=1,col=1,xlab='d',ylab='')
lines(d,sens.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,sens.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
grid()

## Scenario 3
load('sim_scen3_M1000_B5000.RData')
matplot(ascissa,t(y_i),type='l',ylab='',col=gray((1-x_i)*0.6+0.2),lwd=1,lty=1,main='',ylim=c(-2,6.5),xlab='')
title(ylab = "Scenario 3", mgp = c(1.5, 0.4, 0))

plot(d,FWER.unadj,ylim=c(0,1),main='',type='l',lwd=2,lty=1,pch=1,col=1,
     ylab='',xlab='d')
lines(d,FWER.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,FWER.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
grid()
abline(h=0.05,lty=1)

plot(d,FDR.unadj,ylim=c(0,1),main='',type='l',lwd=2,lty=1,pch=1,col=1,xlab='d',ylab='')
lines(d,FDR.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,FDR.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
grid()
abline(h=0.05,lty=1)

plot(d,spec.unadj,ylim=c(0,1),main='',type='l',lwd=2,lty=1,pch=1,col=1,xlab='d',ylab='')
lines(d,spec.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,spec.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
grid()

plot(d,sens.unadj,ylim=c(0,1),main='',type='l',lwd=2,lty=1,pch=1,col=1,xlab='d',ylab='')
lines(d,sens.Fmax,col='gray80',lwd=2,lty=1,pch=1,type='l')
lines(d,sens.FDR,col='gray60',lwd=2,lty=1,pch=1,type='l')
grid()


par(mar=c(0.2,4,1,0.5))
plot(-1:1,-1:1,col=0,ylab='',xlab='',xaxt='n',yaxt='n',bty = 'n')
#legend('center',legend=c('unadj','Fmax'),text.col=c(col.blu[3],col.rosso[3]),horiz=TRUE,text.font=2,cex=2,bty = 'n')
legend('center',legend=c('unadj','fBH','Fmax'),
       col=c('black','gray60','gray80'),
       text.font=2,cex=1.2,lty=1,lwd=2,horiz=TRUE)


FDR.FDR - 1.96*sqrt(FDR.FDR*(1-FDR.FDR)/1000)
