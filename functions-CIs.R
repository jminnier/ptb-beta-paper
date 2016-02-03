#Authors: Jessica Minnier, Lu Tian, Tianxi Cai
#
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

CI.Density <- function(bptb,tsd=NULL,alp=0.95)
{
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


#takes whole matrix of rbind(btrue, bhat, sigmahat, tmpptb)
simulci <- function(bptb,type="Resample",sighat=0,nn=100,rtn="other") {
  
  btrue = bptb[1,]
  bhat = bptb[2,]
  bsig = bptb[3,]
  bptb = bptb[-(1:3),]
  stand = abs(bptb - VTM(bhat,nrow(bptb)))/VTM(bsig+sighat/nn,nrow(bptb))
  #tp0 = apply(bptb,2,function (k) sum(k==0)/length(k))
  tp0 = apply(bptb==0,2,mean)
  if (type=="OLS") {ci.ind = NULL; stand.max.j = apply(stand,1,max)}else{
    ci.ind = which(tp0>=p0.hi)
    if(length(stand[,tp0<p0.hi,drop=F])>0) {stand.max.j = apply(stand[,tp0<p0.hi,drop=F],1,max)}else{
      stand.max.j = 0
    }
  }
  gamma = quantile(stand.max.j,.95)
  
  bcenter = bhat
  bcenter[ci.ind] = 0
  bsig2 = bsig
  bsig2[ci.ind] = 0
  tmpCI = rbind((bcenter - gamma*bsig2),(bcenter + gamma*bsig2))
  tmpCovp = prod((tmpCI[1,]<=btrue)*(tmpCI[2,]>=btrue))
  tmpLen = sum(tmpCI[2,] - tmpCI[1,])
  
  
  tmpout = c(tmpCovp, tmpLen, gamma)
  
  if(rtn=="CI") {tmpout = tmpCI}
  
  tmpout
}


#tmpCImat.simul = simulci.olsgamma(beta.true,tmpbhat,tmpbhat.ini,tmpsd,tmpsd.ols,tmpptb[,-1],iniptb)
simulci.olsgamma <- function(btrue,bhat,bhat.ini,bsig,bsig.ini,bptb,bptb.ini,type="Resample") {
  
  stand = abs(bptb - VTM(bhat,nrow(bptb)))/VTM(bsig,nrow(bptb))
  stand.ini = abs(bptb.ini - VTM(bhat.ini,nrow(bptb.ini)))/VTM(bsig.ini,nrow(bptb.ini))
  tp0 = apply(bptb==0,2,mean)
  if (type=="OLS") {ci.ind = NULL; stand.max.j = apply(stand.ini,1,max)}else{
    ci.ind = which(tp0>=p0.hi)
    if(length(stand[,tp0<p0.hi,drop=F])>0) {stand.max.j = apply(stand.ini[,tp0<p0.hi,drop=F],1,max)}else{
      stand.max.j = 0
    }
  }
  gamma = quantile(stand.max.j,.95)
  
  bcenter = bhat
  bcenter[ci.ind] = 0
  bsig2 = bsig
  bsig2[ci.ind] = 0
  tmpCI = rbind((bcenter - gamma*bsig2),(bcenter + gamma*bsig2))
  tmpCovp = prod((tmpCI[1,]<=btrue)*(tmpCI[2,]>=btrue))
  tmpLen = sum(tmpCI[2,] - tmpCI[1,])
  
  
  c(tmpCovp, tmpLen, gamma)
}	

simulci.alp <- function(bhat,bsig,bptb) {
  stand = abs(bptb - VTM(bhat,nrow(bptb)))/VTM(bsig,nrow(bptb))
  tp0 = apply(bptb==0,2,mean)
  ci.ind = which(tp0>=p0.hi)
  if(length(stand[,tp0<p0.hi,drop=F])>0) {stand.max.j = apply(stand[,tp0<p0.hi,drop=F],1,max)}else{
    stand.max.j = 0
  }
  gamma = quantile(stand.max.j,.95)
  alp = 2*(1-pnorm(gamma))
  alp
}

#Highest Density Region CI function
#apply this function of a matrix of perturbed samples with a vector of bias and true betas at the top     

myhdrfun.ci.new <- function(bptb, biasshift=FALSE,alp=0.05,rtn="other") {
  failurecount.ci = 0; 	failurecount.ci.NA = 0
  bias = bptb[1]; btrue = bptb[2];	bptb = bptb[-(1:2)]
  tp0 = sum(bptb==0)/length(bptb);	tptb = bptb[bptb!=0]
  tmpCI = 1
  ##if(tp0 < p0.no) {tmpCI = c(0, hdr(tptb,prob=95,h=bw.nrd(tptb))$hdr)}else{
  if(tp0 > p0.hi){tmpCI = c(1,0,0)}else{
    if (tp0 < max(p0.lo,alp)) {
      myprob = 100*(1-alp+alp*(tp0+p0.lo)*(tp0 > alp))
      tmpHDR = hdr(tptb,h=bw.nrd(tptb),prob=myprob)
      tmpCI = c(1*(tp0 > alp), tmpHDR$hdr+ biasshift*bias)
    }else{
      myprob = ((1-alp-tp0)*100)/(1-tp0)
      tmpHDR = hdr(tptb,h=bw.nrd(tptb),prob=myprob)
      tmpCI = c(1, tmpHDR$hdr + biasshift*bias)
    }
  }
  
  #print(tmpCI)		
  tmpCI.HDR = tmpCI
  tmpCovp0 = tmpCI[1];  tmpCI = tmpCI[-1]
  ii = NULL
  if(length(tmpCI)>0) {ii = seq(1, length(tmpCI), by=2)}
  
  tmpCovp = sum((tmpCI[ii]<=btrue)*(tmpCI[ii+1]>=btrue)) + tmpCovp0*(btrue==0)
  tmpCovp0 = tmpCovp0 + sum((tmpCI[ii]<=0)*(tmpCI[ii+1]>=0))
  tmpLen = sum(tmpCI[ii+1] - tmpCI[ii])
  
  tmpCovp = 1*(tmpCovp > 0);	tmpCovp0 = 1*(tmpCovp0==0)
  output = c(tmpCovp, tmpCovp0, tmpLen, failurecount.ci, failurecount.ci.NA)
  if(rtn=="CI") {output = c(failurecount.ci, failurecount.ci.NA, tmpLen,tmpCovp0, tmpCI.HDR)}
  output
  ##c(tmpCovp, tmpCovp0, tmpLen, tmpCI)
}



myhdrfun.ci.new2 <- function(bptb, biasshift=FALSE,alp=0.05,rtn="other") {
  failurecount.ci = 0; 	failurecount.ci.NA = 0
  bias = bptb[1]; btrue = bptb[2];	bptb = bptb[-(1:2)]
  tp0 = sum(bptb==0)/length(bptb);	tptb = bptb[bptb!=0]
  tmpCI = 1
  ##if(tp0 < p0.no) {tmpCI = c(0, hdr(tptb,prob=95,h=bw.nrd(tptb))$hdr)}else{
  if(tp0 > p0.hi){tmpCI = c(1,0,0)}else{
    if (tp0 < max(p0.lo,alp)) {
      myprob = 100*(1-alp+alp*(tp0+p0.lo)*(tp0 > alp))
      tmpHDR = hdr(tptb,h=bw.nrd(tptb),prob=myprob)
      tmpCI = c(1*(tp0 > alp), tmpHDR$hdr+ biasshift*bias)
    }else{
      myprob = ((1-alp-tp0)*100)/(1-tp0)
      tmpHDR = hdr(tptb,h=bw.nrd(tptb),prob=myprob)
      tmpCI = c(1, tmpHDR$hdr + biasshift*bias)
    }
  }
  
  #print(tmpCI)		
  tmpCI.HDR = tmpCI
  tmpCovp0 = tmpCI[1];  tmpCI = tmpCI[-1]
  ii = NULL
  if(length(tmpCI)>0) {ii = seq(1, length(tmpCI), by=2)}
  
  tmpCovp = sum((tmpCI[ii]<=btrue)*(tmpCI[ii+1]>=btrue)) + tmpCovp0*(btrue==0)
  tmpCovp0 = tmpCovp0 + sum((tmpCI[ii]<=0)*(tmpCI[ii+1]>=0))
  tmpLen = sum(tmpCI[ii+1] - tmpCI[ii])
  
  tmpCovp = 1*(tmpCovp > 0);	tmpCovp0 = 1*(tmpCovp0==0)
  output = c(tmpCovp, tmpCovp0, tmpLen, failurecount.ci, failurecount.ci.NA)
  if(rtn=="CI") {output = c(failurecount.ci, failurecount.ci.NA, tmpLen,tmpCovp0, tmpCI.HDR)}
  output
  ##c(tmpCovp, tmpCovp0, tmpLen, tmpCI)
}

#Dec 17, 2010, parallel Normal CI as HDR for \hat{c_2}
mynormfun.ci <- function(bptb) {
  bhat = bptb[1]; bsd=bptb[2]; bptb= bptb[-(1:2)]
  tp0 = mean(bptb==0)
  if(tp0 >p0.hi) {tmpCI = c(0,0)}else{
    if(tp0 < max(p0.lo,.05)) {
      myprob = .95+.05*(tp0+p0.lo)*(tp0 > 0.05)
      print(myprob)
      myz = -1*qnorm((1-myprob)/2)
      tmpCI = c(bhat - myz*bsd, bhat+myz*bsd)
    }else{
      tmpCI = c(bhat - 1.96*bsd, bhat+1.96*bsd)
    }
  }
  tmpCI
}	


#counts how many of the true signals are included in the model
truemodel = function(bhat, btrue) {
  ind.non0 = which(btrue!=0)
  sum(bhat[ind.non0]!=0)
}
