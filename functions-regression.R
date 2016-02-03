#Author: Jessica Minnier, Lu Tian, Tianxi Cai

##============================================================================##
#     Code to implement simulations and data analyses from 
# "A perturbation method for inference on regularized regression estimates"
# Copyright (C) 2010-2016  Jessica Minnier <minnier-at-ohsu.edu>
#   
#   This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
##============================================================================##


library(MASS)
library(lars)
library(lasso2) 
library(hdrcde)

VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}


SIM.FUN = function(nn,beta0=beta.true,rho=0,sig0=mysig0)
{
  p0=length(beta0)
  Xi=mvrnorm(nn,rep(0,p0),Sigma=diag(1-rho,p0)+rho);
  Yi=rnorm(nn)*sig0 + Xi%*%beta0
  cbind(Yi,Xi)
}

Dhat.LASSO = function(data, K0.fold=1,s0=NULL)
{
  if(K0.fold>1)
  {
    Dhat = 0; nn = nrow(data); nv = floor(nn/K0.fold)
    for(k in 1:K0.fold)
    {
      indv = 1:nv + (k-1)*nv; indt = setdiff(1:nn,indv)
      beta.t = Est.LASSO(data[indt,])$beta
      Dhat = Dhat + mean(abs(data[indv,1]-cbind(1,data[indv,-1])%*%beta.t))
    }
    Dhat = Dhat/K0.fold
  }else{
    betahat = Est.LASSO(data,s0=s0)$beta
    Dhat = mean(abs(data[,1]-cbind(1,data[,-1])%*%betahat))	  }
  Dhat  	
} 

Est.LASSO = function(data,s0=NULL,type0=NULL,rtn="ALL")
{
  y = data[,1]; x = data[,-1]
  tmpfit = lars(x,y,type="lasso",normalize = F, intercept=TRUE)
  if(is.null(s0))
  {
    #set.seed(315)
    fit.cv  = cv.lars(x,y,K=5,plot.it=F,normalize = F, intercept=TRUE, se=F)
    s0  = min(fit.cv$fraction[fit.cv$cv==min(fit.cv$cv)])
    type0 = "fraction"
  }
  betahat  = coef(tmpfit,s=s0,mode=type0); #print(c(s0,betahat))
  alphahat = mean(y) - sum(apply(x,2,mean)*betahat)
  list("beta"=c(alphahat,betahat),"s0"=s0,"type0"=type0)
}

Est.ALASSO1 = function(data,s0=NULL,type0=NULL,rtn="ALL")
{
  y = data[,1]; x = data[,-1]
  w = 1/abs(lm(y~x-1)$coef); x.t = x/VTM(w,nrow(x))
  tmpfit = lars(x.t,y,type="lasso",normalize = F, intercept=F)
  if(is.null(s0))
  {
    fit.cv  = cv.lars(x.t,y,K=5,plot.it=F,normalize = F, intercept=F)
    s0  = fit.cv$fraction[fit.cv$cv==min(fit.cv$cv)]
    type0 = "fraction"
  }
  betahat  = coef(tmpfit,s=s0,mode=type0)/w
  list("beta"=betahat,"s0"=s0,"type0"=type0)
}

Est.ALASSO2 = function(data,s0=NULL,rtn="ALL")
{
  y = data[,1]; x = data[,-1]; nn=length(y)
  beta.ini = lm(y~x)$coef; w = 1/abs(beta.ini[-1]); x.t = x/VTM(w,nrow(x))
  m0 = 500; lam.all = s0.all = (1:m0)/m0
  fit.all = l1ce(y~x.t, standardize=F, bound=s0.all)
  beta.all = coef(fit.all)/VTM(c(1,w),m0)
  df.all = apply(beta.all[,-1]!=0,1,sum)
  Ahat = t(cbind(1,x))%*%cbind(1,x)/nn 
  BIC.lam = diag((beta.all-VTM(beta.ini,m0))%*%Ahat%*%t(beta.all-VTM(beta.ini,m0))) + log(nn)*df.all/nn
  m.opt = (1:m0)[BIC.lam==min(BIC.lam)]
  fit.hat = fit.all[[m.opt]]; fit.hat$coef = fit.hat$coef/c(1,w)
  return(c("beta"=beta.all[m.opt,], "lam0" =fit.hat$Lagrangian, "b0"=fit.hat$bound))
}

L2Norm <- function(data,coef,intercept=F)
{
  yy = data[,1]; xx = data[,-1]; if(intercept){xx=cbind(1,xx)}
  apply((yy - xx%*%t(coef))^2,2,mean)
}

Est.ALASSO = function(data,s0=NULL,rtn="ALL",intercept=F,s.opt=NULL, myratio=1,bictype="min")
{
  y = data[,1]; x = data[,-1]; nn=length(y); Ahat = t(x)%*%x/nn; pp=ncol(data)
  fit.ini = lm(y~x-1)
  beta.ini = fit.ini$coef; #use beta.hat.lambda 
  #sighat2 = (nn/(nn-pp))*(summary(fit.ini)$sigma)^2
  sighat2 = (summary(fit.ini)$sigma)^2
  #beta.ini = lm.ridge(y~x-1,lambda=min(nn^.01,log(nn)))$coef
  w = 1/abs(beta.ini);  x.w = x/VTM(w,nrow(x))
  tmpfit = lars(x.w,y,type="lasso",normalize = F, intercept=F)
  
  ##*********************************************************##
  ## nlminb is meant for functions that are smooth and might ##
  ## give local min for non-smooth functions like BIC.lam    ##       
  ##*********************************************************##
  
  lam.all = range(log(tmpfit$lambda+0.0001))
  lam.all = exp(seq(lam.all[1],lam.all[2],length=2000))
  coef.all = t(coef(tmpfit,s=lam.all,mode="lambda"))/w  # pz x n.lam 
  BIC.lam = diag(t(coef.all-beta.ini)%*%Ahat%*%(coef.all-beta.ini)) +
    lam.all*apply( (coef.all-beta.ini)^2*(coef.all!=0)/(abs(beta.ini)*pmax(abs(coef.all),abs(beta.ini)/sqrt(nn))),2,sum)/nn
  bicn = min(log(nn),nn^0.1)
  if(bictype=="logn") bicn=log(nn)
  BIC.lam = BIC.lam + sighat2*bicn/nn*apply(coef.all!=0,2,sum)
  ind.opt  = 1 + which.min(BIC.lam[-1]);   lam.opt = lam.all[ind.opt]; 
  
  betahat = coef(tmpfit,s=lam.opt,mode="lambda")/w
  beta.out = c("lam0" = lam.opt, "beta"=betahat*myratio,  "bols"=beta.ini*myratio)
  if(rtn=="ALL")
  { 
    biashat = lam.opt*solve(Ahat)%*%sign(betahat)/abs(betahat)/nn
    biashat[abs(biashat)==Inf]=0
    beta.out = c(beta.out,biashat); 
    return(list("beta.out"=beta.out,"Ahat"=Ahat,"sighat2"=sighat2))
  }else{
    return(c("beta.out"=beta.out)) }
}



CV.FUN <- function(data)
{
  y = data[,1]; x = data[,-1]; nn=length(y)
  w = 1/abs(lm(y~x)$coef[-1]); x.w = x/VTM(w,nrow(x))
  s0.all = (1:100)/100
  fit.all = l1ce(y~x.w, standardize=F, bound=s0.all)
  K = 5; L2cv = NULL
  for(k in 1:K)
  {
    indv = 1:floor(nn/K) + (k-1)*floor(nn/K)
    indt = setdiff(1:nn,indv)
    fitk = l1ce(y~x.w,subset=indt, standardize=F, bound=s0.all)
    L2cv = rbind(L2cv, L2Norm(cbind(y,x.w)[indv,],coef(fitk),intercept=T))
  }
  L2cv = apply(L2cv,2,mean)
  s0 = min(s0.all[L2cv==min(L2cv)])
  betahat  = l1ce(y~x.w, standardize=F, bound=s0)
  list("beta"=coef(betahat)/w, "s0"=s0, "lam0" =betahat$Lagrangian,"b0"=betahat$bound)    
}

###############################################################################################################################
###############################################################################################################################
##################### Copied functions from Summary-LASSO-jm-new                   ############################################
##################### February 9, 2010                                             ############################################
###############################################################################################################################
###############################################################################################################################



CI.Density <- function(bptb,tsd=NULL,alp=0.95)
{
  #browser()
  tp0=apply(bptb==0,2,mean); tptb=bptb; tptb[bptb==0]=NA
  tmu=apply(tptb,2,mean,na.rm=T); if(is.null(tsd)){tsd=apply(tptb,2,sd,na.rm=T)}
  tsd[is.na(tsd)]=0
  tmu[is.na(tmu)]=0
  CovP <- function(cc,mu,sd,p0)
  {
    const=-2*log(cc/(1-p0)*sqrt(2*pi)*sd)
    if(p0<1) {p0*(p0>=cc)+(1-p0)*pchisq(const,df=1)-alp}
    else p0*(p0>=cc)
  }
  out=tc0=NULL
  for(ip in 1:length(tp0))
  {
    tc0=c(tc0,uniroot(CovP,interval=c(1e-5,1000),mu=tmu[ip],sd=tsd[ip],p0=tp0[ip])$root)
    tconst = sqrt(-2*log(tc0[ip]/(1-tp0[ip])*sqrt(2*pi)*tsd[ip]))
    tmpout = ifelse(tp0[ip]>=tc0[ip],0,NA)
    tmpout = c(tmpout, tmu[ip]-tconst*tsd[ip], tmu[ip]+tconst*tsd[ip])
    out = rbind(out,tmpout)
  }
  out
}


naivesd <- function(b.hat, A.hat, n) {
  not0 = which(b.hat!=0)
  sd = rep(0,length(b.hat))
  if(length(not0)>0) {
    A.hat = matrix(A.hat, ncol=length(b.hat));		A.hat11 = A.hat[not0,not0]
    sd[not0] = sqrt(diag(solve(A.hat11))/n)
  }	
  sd
}

#took out optimization of quantile length to run faster
CI.quantlen.vec <- function(bvec,alp=seq(0,0.05,by=1/500)) {
  c(1,2,3, quantile(bvec, c(0.975,0.025)))
}

CI.quantlen.vec.OLD <- function(bvec,alp=seq(0,0.05,by=1/500)) {
  CI.l = quantile(bvec, alp)
  CI.h = quantile(bvec, alp+.95)
  alp.ind = which.min(CI.h-CI.l)
  c(alp[alp.ind], CI.h[alp.ind], CI.l[alp.ind], quantile(bvec, c(0.975,0.025)))
}

CI.quantlen <- function(bptb,bhat=NULL,alp=seq(0,0.05,by=1/500)) {
  if(length(bhat)==0) {
    output = t(apply(bptb, 2, CI.quantlen.vec))
    cbind(output[,1], output[,3], output[,2], output[,5], output[,4])
  }
  else {
    bvec = bptb - VTM(bhat,length(bptb[,1])) #subtract row by row
    output = t(apply(bvec, 2, CI.quantlen.vec))
    cbind(output[,1], bhat - output[,2:5])
  }
}         

myhdrfun <- function(colvec) {
  if(sd(colvec[-1])==0) {0}
  else {hdr(colvec[-1],h=colvec[1])$mode}
}

myhdrfun.mode.FAKE <- function(bptb) {0}		
myhdrfun.mode <- function(bptb) {
  failcount = 0
  num0 = length(bptb) - sum(bptb==0) 
  hdrmode = 0
  if(num0>3) {
    tp0 = sum(bptb==0)/length(bptb)
    if(tp0<.5) {
      tptb = bptb[bptb!=0]
      
      tdens = try(density(tptb, bw=hdrbw(tptb,.1)),silent=T)
      if(class(tdens)=="try-error") {tdens = density(tptb); failcount = 1}
      cmax = max(tdens$y)*(1-tp0)
      
      hdrmode = (cmax>tp0/(1-tp0))*hdr(tptb,h=tdens$bw)$mode
    }
  }
  
  c(hdrmode, failcount) #failcount is 1 if the hdrbw fails
}

myhdrfun.ci.FAKE <- function(bptb) { c(1,2,3) }

#apply to matrix of beta.true, bhat.densmode, tmpptb
myhdrfun.ci.WITHMODE <- function(bptb) {
  #col of c(btrue, bptb)
  btrue = bptb[1]
  bmode = bptb[2]
  bptb = bptb[-(1:2)]
  
  tp0 = sum(bptb==0)/length(bptb)
  tptb = bptb[bptb!=0]
  
  
  if (bmode == 0) {
    tmpCI = 1
    if (tp0 < 0.95) {
      
      if(tp0 > 0.5) {
        bw.opt = density(tptb)$bw
        tmpCI = c(1, hdr(tptb,h=bw.opt,100*(.95-tp0)/(1-tp0))$hdr) #divide by 1-tp0???
      }
      else {
        prob.seq = sort(seq(1,99,length=500),decreasing=T)
        covp=95
        ii=0
        while(covp>0) {
          ii = ii + 1
          tmpHDR = hdr(tptb,prob=prob.seq[ii])
          cc = tmpHDR$falpha
          covp = tp0*((tp0/(1-tp0))>=cc) + (prob.seq[ii]/100)*(1-tp0) - 0.95
        }
        tmpCI = c(1*(tp0/(1-tp0)>=cc), tmpHDR$hdr)
      }
    }
  }
  else { tmpCI = c(0, hdr(bptb,prob=.95)$hdr) }	
  
  
  tmpCovp0 = tmpCI[1]
  tmpCI = tmpCI[-1]
  
  ii = NULL
  if(length(tmpCI)>0) {ii = seq(1, length(tmpCI), by=2)}
  
  tmpCovp = sum((tmpCI[ii]<=btrue)*(tmpCI[ii+1]>=btrue)) + tmpCovp0*(btrue==0)
  tmpCovp0 = tmpCovp0 + sum((tmpCI[ii]<=0)*(tmpCI[ii+1]>=0))
  tmpLen = sum(tmpCI[ii+1] - tmpCI[ii])
  
  tmpCovp = 1*(tmpCovp > 0)
  tmpCovp0 = 1*(tmpCovp0==0)
  c(tmpCovp, tmpCovp0, tmpLen)
  #Covp, Cov0, Len 
  
}

#apply to matrix of beta.true, tmpptb
myhdrfun.ci <- function(bptb,rtn="coverage") {
  failurecount.ci = 0
  failurecount.ci.NA = 0
  
  #col of c(btrue, bptb)
  btrue = bptb[1]
  bptb = bptb[-1]
  
  tp0 = sum(bptb==0)/length(bptb)
  tptb = bptb[bptb!=0]
  
  tmpCI = 1
  if (tp0 < 0.95) {
    
    if(tp0 > 0.5) {
      bw.opt = density(tptb)$bw
      tmpCI = c(1, hdr(tptb,h=bw.opt,100*(.95-tp0)/(1-tp0))$hdr) #divide by 1-tp0???
    }
    else {
      prob.seq = sort(seq(1,99,length=500),decreasing=T)
      covp=95
      ii=0
      while(covp>0) {
        ii = ii + 1
        tmpHDR = try(hdr(tptb,prob=prob.seq[ii]),silent=T)
        if(class(tmpHDR)=="try-error") {tmpHDR = hdr(tptb,h=density(tptb)$bw,prob=prob.seq[ii]); failurecount.ci = 1}
        cc = tmpHDR$falpha
        covp = tp0*((tp0/(1-tp0))>=cc) + (prob.seq[ii]/100)*(1-tp0) - 0.95
      }
      tmpCI = c(1*(tp0/(1-tp0)>=cc), tmpHDR$hdr)
      if (length(tmpHDR$hdr)<2) {tmpCI = c(1*(tp0/(1-tp0)>=cc), 0, 0); failurecount.ci.NA = 1}
    }
  }
  
  tmpCI.HDR = tmpCI	
  tmpCovp0 = tmpCI[1]
  tmpCI = tmpCI[-1]
  
  ii = NULL
  if(length(tmpCI)>0) {ii = seq(1, length(tmpCI), by=2)}
  
  tmpCovp = sum((tmpCI[ii]<=btrue)*(tmpCI[ii+1]>=btrue)) + tmpCovp0*(btrue==0)
  tmpCovp0 = tmpCovp0 + sum((tmpCI[ii]<=0)*(tmpCI[ii+1]>=0))
  tmpLen = sum(tmpCI[ii+1] - tmpCI[ii])
  
  tmpCovp = 1*(tmpCovp > 0)
  tmpCovp0 = 1*(tmpCovp0==0)
  output = c(tmpCovp, tmpCovp0, tmpLen, failurecount.ci, failurecount.ci.NA)
  #Covp, Cov0, Len 
  if(rtn=="CI") {output = c(failurecount.ci, failurecount.ci.NA, tmpLen,tmpCovp0, tmpCI.HDR)}
  output	
}

#counts how many of the true signals are included in the model
truemodel = function(bhat, btrue) {
  ind.non0 = which(btrue!=0)
  sum(bhat[ind.non0]!=0)
}
