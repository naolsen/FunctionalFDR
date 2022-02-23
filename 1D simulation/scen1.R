

library(fda)

rmvunif = function(n,p,min=0,max=1){
  tmp = as.matrix(rep(n,p))
  data = apply(tmp,1,runif,min=min,max=max)
  return(data)
}
rmvexp = function(n,p,lambda=1){
  tmp = as.matrix(rep(n,p))
  data = apply(tmp,1,rexp,rate=lambda) - 1
  return(data)
}
rmvgauss = function(n,p,mu=0,sigma=1){
  tmp = as.matrix(rep(n,p))
  data = apply(tmp,1,rnorm,mean=mu,sd=sigma)
  return(data)
}


# Code for Fmax function:
source('Fmax.R')

# Simulation --------------------------------------------------------------

# Continuous covariate ----------------------------------------------------
scen = 1

d = seq(0,5,len=6) # len = 11
M = 1000
B = 5000
beta0 = 0
n = 10 
ncoef = 40
J = 200
ascissa = seq(0,1,len=J)
Incoef = diag(ncoef)
basis = create.bspline.basis(range(ascissa),nbasis=ncoef,norder=4)

mu = 0
sigma = 1
lambda = 1
sigmacoeff_indip = sigma^2 * diag(ncoef)


if(scen==1){
  h=ncoef/4
}else if(scen==2){
  h=ncoef/2
}else if(scen==3){
  h=ncoef*3/4
}
mucoef = c(rep(1,h),rep(0,ncoef-h))
fd_funct = fd(coef=mucoef,basisobj=basis) 
funct = as.vector(eval.fd(ascissa,fd_funct))

x_i = seq(0,1,len=n)


spec.unadj = spec.Fmax = spec.FDR = numeric(length(d))
FWER.unadj = FWER.Fmax = FWER.FDR = numeric(length(d))
FDR.unadj = FDR.Fmax = FDR.FDR = numeric(length(d))
sens.unadj = sens.Fmax = sens.FDR = numeric(length(d))
pval.unadj = pval.Fmax = pval.FDR = array(dim=c(M,J,length(d)))

alpha = 0.05



for(index in 1:length(d)){
  set.seed(15072017)
  print(index)
  d.tmp = d[index]
  beta1 = d.tmp* funct
  null = which(beta1 == 0)
  alt = which(beta1 != 0)
  for(sim in 1:M){
    c_e_i = rmvgauss(n,ncoef,mu,sigma)
    f_e_i = fd(coef=t(c_e_i),basisobj=basis) 
    e_i = t(eval.fd(ascissa,f_e_i))
    # data
    y_i = matrix(nrow=n,ncol=J)
    for(ii in 1:n){
      y_i[ii,] = beta0 + beta1*x_i[ii] + e_i[ii,] 
    }
    Fmax_res = Fmax(y_i ~ x_i,B=B,method='responses')
    pval.Fmax[sim,,index] = Fmax_res$adjusted_pval_F
    pval.unadj[sim,,index] = Fmax_res$unadjusted_pval_F
    pval.FDR[sim,,index] = p.adjust(pval.unadj[sim,,index],method='BH')
  }
  if(length(alt)>0){
    sens.unadj[index] = mean(rowSums(pval.unadj[,alt,index] < alpha)/length(alt))
    sens.Fmax[index]  = mean(rowSums(pval.Fmax[,alt,index] < alpha)/length(alt))
    sens.FDR[index]   = mean(rowSums(pval.FDR[,alt,index] < alpha)/length(alt))
    
  }
  if(length(null)>0){
    spec.unadj[index]  = mean(rowSums(pval.unadj[,null,index] < alpha)/length(null))
    spec.Fmax[index] = mean(rowSums(pval.Fmax[,null,index] < alpha)/length(null))
    spec.FDR[index]  = mean(rowSums(pval.FDR[,null,index] < alpha)/length(null))
    
    FWER.unadj[index]  = mean((rowSums(pval.unadj[,null,index] < alpha))>0)
    FWER.Fmax[index] = mean((rowSums(pval.Fmax[,null,index] < alpha))>0)
    FWER.FDR[index]  = mean((rowSums(pval.FDR[,null,index] < alpha))>0)
    
    FDP.unadj = rowSums(pval.unadj[,null,index] < alpha)/rowSums(pval.unadj[,,index] < alpha)
    FDP.unadj[which(is.nan(FDP.unadj))] = 0
    FDP.Fmax  = rowSums(pval.Fmax[,null,index] < alpha)/rowSums(pval.Fmax[,,index] < alpha)
    FDP.Fmax[which(is.nan(FDP.Fmax))] = 0
    FDP.FDR   = rowSums(pval.FDR[,null,index] < alpha)/rowSums(pval.FDR[,,index] < alpha)
    FDP.FDR[which(is.nan(FDP.FDR))] = 0
    
    FDR.unadj[index]   = mean(FDP.unadj)
    FDR.Fmax[index]  = mean(FDP.Fmax)
    FDR.FDR[index]   = mean(FDP.FDR)
  }
  
}

save.image(paste0('sim_scen',scen,'_M',M,'_B',B,'.RData'))




