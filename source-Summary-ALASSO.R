#Author: Jessica Minnier

## ********************************************************************************** ## 
## **** File takes alasso simulation output and creates summary csv file         **** ##                
## ********************************************************************************** ## 
## ********************************************************************************** ##
## 
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

replaceSum = FALSE

source('functions-CIs.R')

################################################################################################################
############## Ahat is Ahat + diag P(betaptb!=0)/ols est 
################################################################################################################


junklam  = junk[,1]; 
junkhat  = junk[,1+1:pz]; 				#alasso estimates
junkini  = junk[,1+pz+1:pz];			#ols estimates
junkbias = junk[,1+2*pz+1:pz]; 			#bias estimate
junkAhat = junk[,1+3*pz+1:(pz^2)]		#X'X/n
junksig2 = junk[,1+3*pz+(pz^2)+1]		#sighat2
junkptb  = junk[,-(1:(2+3*pz+pz^2))]	#ptbed

beta.mean <- apply(junkhat,2,mean)
vec0 <- rep(0,pz)



(nptb0 <- length(junk[,1]))

tmpcsv = paste(outputwd.sum,"summary.b",sprintf("%1d",beta.true[4]*10),"n",sampsize,"ptb",nptb0,"pz",pz,"sig",mysig0,"corX",cor.X,".BIC.csv",sep="")
if((!file.exists(tmpcsv))||replaceSum){ #if file doesn't exist or replaceSum =TRUE
  
  
  
  
  #plot(density(junklam),main="Lambda Density, log(n)*n^.1")
  
  junkbias2=junkp0= junksd0 = junksd =junksd2= junksd.ols=
    junkcov.ols = junkcov.naive = junkcov.norm = junkcov.normbc = junkcov.dens = junkcov.densbc = 
    junkcov.ols.mean = junkcov.naive.mean = junkcov.norm.mean = junkcov.normbc.mean = junkcov.dens.mean = junkcov.densbc.mean =
    junklen.ols = junklen.naive = junklen.norm = junklen.dens = junksd.naive =
    junkcov.norm2 = junkcov.dens2 = junkcov.normbc2 = junkcov.densbc2 = junkcov.norm.mean2 = junkcov.dens.mean2 = 
    junkcov.densbc.mean2 = junkcov.normbc.mean2 = junklen.norm2 = junklen.dens2 = 
    junkcov.dens.null = junkcov.dens.mean.null = junkcov.densbc.null = junkcov.densbc.mean.null = junklen.dens.null = junkmodel=
    junkcov0.ols = junkcov0.naive= junkcov0.norm2 = junkcov0.norm = junkcov0.normbc2 = junkcov0.normbc = junkcov0.dens2 = junkcov0.dens=
    junkcov0.dens.null = junkcov0.densbc2 = junkcov0.densbc = junkcov0.densbc.null =junkcov.quantlen = junkcov.95 = junkcov.quantlen.mean =
    junkcov.95.mean = junkcov.quantlen0 = junkcov.950 = junklen.95 = junklen.quantlen = junkAlp.opt = 
    junkcov.quantlenBC = junkcov.95BC = junkcov.quantlen0BC = junkcov.950BC = 
    junkcov.quantlen.star = junkcov.95.star = junkcov.quantlen.mean.star =
    junkcov.95.mean.star = junkcov.quantlen0.star = junkcov.950.star = junklen.95.star = junklen.quantlen.star = junkAlp.opt.star = 
    junkcov.quantlenBC.star = junkcov.95BC.star = junkcov.quantlen0BC.star = junkcov.950BC.star = 
    junklen.normA = junkcov.normA = junkcov.normA.mean = junkcov0.normA = 
    junkcov.normbcA = junkcov.normbcA.mean = junkcov0.normbcA = 
    junklen.normA.ZouAhat = junkcov.normA.ZouAhat = junkcov.normA.mean.ZouAhat = junkcov0.normA.ZouAhat = 
    junkcov.normbcA.ZouAhat = junkcov.normbcA.mean.ZouAhat = junkcov0.normbcA.ZouAhat = 
    junkbias2.ZouAhat= junksd2.ZouAhat = junksd.ZouAhat = 
    junkbhat.densmode = 
    junklen.norm.densmode =  junkcov.norm.densmode = junkcov.norm.mean.densmode =  junkcov0.norm.densmode=
    junkmodel.densmode= junkp0.densmode = junkcov.hdr =  junkcov0.hdr = junklen.hdr = 
    junkcov.hdr0 =  junkcov0.hdr0 = junklen.hdr0 =
    junkcov.hdr1 =  junkcov0.hdr1 = junklen.hdr1 =
    junkcov.hdr2 =  junkcov0.hdr2 = junklen.hdr2 =
    junkcov.hdr3 =  junkcov0.hdr3 = junklen.hdr3 =
    junkcov.hdr4 =  junkcov0.hdr4 = junklen.hdr4 =
    junkbhat.densmedian = 
    junklen.norm.densmedian =  junkcov.norm.densmedian = junkcov.norm.mean.densmedian =  junkcov0.norm.densmedian=
    junkmodel.densmedian= junkp0.densmedian =
    junktruemodel.densmode = junktruemodel.densmedian = junktruemodel = 
    junkbias2.non0 = junkbias2.ZouAhat.non0 = 
    junkbias2.nop0 = junkbias2.nop0.non0 = junksd2.nop0 = 
    junkcov0.normbc2.nop0 = junkcov.normbc2.nop0 = junkcov0.norm.nop0 = junklen.norm.nop0 = junkcov.norm.nop0 = 
    junkcov0.normbc.nop0 = junkcov.normbc.nop0 = 
    junkcov.simul = junklen.simul = junkcov.simul.ols = junklen.simul.ols = junkgam.simul = junkgam.simul.ols = 
    junkcov.norm2p0 = junkcov0.norm2p0 = junklen.norm2p0 = 
    junklam.ptb = junk.p0lo = junk.p0hi = junkcov.norm2infl = junklen.norm2infl = junkcov0.norm2infl = 
    junkcov.hdrBC = junklen.hdrBC = junkcov0.hdrBC =
    junkcov.simul2  = junklen.simul2 = junkgam.simul2 =junkcov.simul3  = junklen.simul3 = junkgam.simul3 =
    junklen.hdr.simul = junkcov.hdr.simul = junksimul.alp = NULL
  
  failurecount.mode = failurecount.ci = failurecount.ci.NA = failurecount.mode0 = failurecount.ci0 = failurecount.ci.NA0 =
    failurecount.mode1 = failurecount.ci1 = failurecount.ci.NA1 = 0
  
  
  for(kk in 1:nrow(junk)){	
    tmpptb = matrix(as.numeric(junkptb[kk,]),ncol=2*pz+1); 
    iniptb = tmpptb[,-(1:(1+pz))]; tmpptb = tmpptb[,1:(1+pz)]
    junklam.ptb = c(junklam.ptb,mean(tmpptb[,1]))
    tmpbhat = as.numeric(junkhat[kk,]); tmpmodel = length(tmpbhat[tmpbhat!=0])
    tmpbhat.ini = as.numeric(junkini[kk,])
    tmpsighat2 = junksig2[kk]
    
    junkmodel = rbind(junkmodel,tmpmodel)		                                                #model size
    dvec = apply(tmpptb!=0,2,mean)[-1]/tmpbhat.ini^2                                           #Ahat = Ahat + p(b*!=0)/bhat^2
    tmpAhatlam = junkAhat[kk,] + (junklam[kk]/tmpn)*as.vector(diag(dvec))
    tmpAhatlam.nop0 = junkAhat[kk,] + (junklam[kk]/tmpn)*as.vector(diag(1/tmpbhat.ini^2))
    
    not0 = which(tmpbhat!=0)				                                                    #Zou2006 Ahat (not Ahat but Jacobian)
    if(length(not0)==0){
      Ahat.submat = 0; tmpsd.ZouAhat = rep(0,pz)
    }else{
      Ahat.submat = tmpn*matrix(as.numeric(junkAhat[kk,]),nrow=pz)[not0,not0,drop=F]
      mat1 = solve(Ahat.submat + 
                     (junklam[kk]/tmpn)*diag(
                       1/(as.numeric(abs(tmpbhat*tmpbhat.ini)[not0])),nrow=length(not0)))
      junkAhat2 = mysig0*(mat1%*%Ahat.submat%*%mat1)
      tmpsd.ZouAhat = rep(0,pz)
      tmpsd.ZouAhat[not0] = sqrt(diag(junkAhat2))
    }
    junksd.ZouAhat = rbind(junksd.ZouAhat, tmpsd.ZouAhat)
    
    ############ emp sd of ptb betas #############################################################
    tmpsd = sqrt(apply( (tmpptb[,-1]-VTM(tmpbhat,nrow(tmpptb)))^2,2,mean));                #tmpsd: 
    junksd = rbind(junksd,tmpsd)                                                                #sd = (ptb - junkhat)^2 ave
    tmpsd.naive = naivesd(tmpbhat, tmpAhatlam,tmpn)
    junksd.naive = rbind(junksd.naive, tmpsd.naive)
    
    tmpCI = cbind(tmpbhat-1.96*tmpsd.naive,tmpbhat+1.96*tmpsd.naive)                  #NaiveCovp, NaiveCovp Mean, NaiveLen
    tmpcov= (tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true)
    tmpcov.mean= (tmpCI[,1]<=beta.mean)*(tmpCI[,2]>=beta.mean)
    junklen.naive=rbind(junklen.naive,tmpCI[,2]-tmpCI[,1])
    junkcov.naive=rbind(junkcov.naive,tmpcov)
    junkcov.naive.mean=rbind(junkcov.naive.mean,tmpcov.mean)
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.naive=rbind(junkcov0.naive,tmpcov0)
    
    tmpsd.ols = sqrt(apply( (iniptb-VTM(tmpbhat.ini,nrow(iniptb)))^2,2,mean));
    junksd.ols=rbind(junksd.ols,tmpsd.ols)                                             #OLSCovp, OLSCovp Mean, OLSLen
    tmpCI=cbind(tmpbhat.ini-1.96*tmpsd.ols,tmpbhat.ini+1.96*tmpsd.ols)      #empirical sd of iniptb
    junklen.ols=rbind(junklen.ols,tmpCI[,2]-tmpCI[,1])
    junkcov.ols=rbind(junkcov.ols,1*(tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true))
    junkcov.ols.mean=rbind(junkcov.ols.mean,1*(tmpCI[,1]<=beta.mean)*(tmpCI[,2]>=beta.mean))
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.ols=rbind(junkcov0.ols,tmpcov0)
    
    tmpp0 = apply(tmpptb[,-1]==0,2,mean); junkp0 = rbind(junkp0,tmpp0)                          #P(beta=0)
    tmpbias=(-1)^(apply(sign(tmpptb[,-1])==1,2,mean)<apply(sign(tmpptb[,-1])==-1,2,mean))/
      max(abs(quantile(iniptb,c(0.075,0.975))))                                 #apply(abs(iniptb),2,quantile,0.95)
    tmpbias = mean(tmpptb[,1])*tmpbias%*%solve(matrix(tmpAhatlam,ncol=pz))/tmpn
    tmpbias = c(tmpbias*(1-tmpp0)) 				       
    tmpbias.non0 = tmpbias*(tmpbhat!=0)                                         #changed from 1-tmpp0
    junkbias2= rbind(junkbias2,tmpbias)
    junkbias2.non0 = rbind(junkbias2.non0, tmpbias.non0)
    tmpsd2 = sqrt(tmpsd^2+tmpbias^2);  			                                                #tmpsd2 includes bias
    junksd2=rbind(junksd2,tmpsd2)
    
    ########## TRY BIAS WITHOUT 1-HAT(P0) TERMS
    tmpbias.nop0=(-1)^(apply(sign(tmpptb[,-1])==1,2,mean)<apply(sign(tmpptb[,-1])==-1,2,mean))/
      max(abs(quantile(iniptb,c(0.075,0.975))))                                     #apply(abs(iniptb),2,quantile,0.95)
    tmpbias.nop0 = mean(tmpptb[,1])*tmpbias.nop0%*%solve(matrix(tmpAhatlam.nop0,ncol=pz))/tmpn
    tmpbias.nop0 = c(tmpbias.nop0*1) 		
    tmpbias.nop0.non0 = tmpbias.nop0*(tmpbhat!=0)                                         #changed from 1-tmpp0
    junkbias2.nop0= rbind(junkbias2.nop0,tmpbias.nop0)
    junkbias2.nop0.non0 = rbind(junkbias2.nop0.non0, tmpbias.nop0.non0)
    tmpsd2.nop0 = sqrt(tmpsd^2+tmpbias.nop0^2);                                                  #tmpsd2 includes bias
    junksd2.nop0=rbind(junksd2.nop0,tmpsd2.nop0)
    
    ############# Bias using Zou's Cov Matrix
    if(length(not0)==0) {tmpbias.ZouAhat = rep(0,pz)
    }else{
      tmpptb.ZouAhat = tmpptb[,-1]
      tmpptb.ZouAhat = tmpptb.ZouAhat[,not0,drop=F]
      tmpbias.ZouAhat =(-1)^(apply(sign(tmpptb.ZouAhat)==1,2,mean)<apply(sign(tmpptb.ZouAhat)==-1,2,mean))/
        max(abs(quantile(iniptb[,not0],c(0.075,0.975))))                                     #apply(abs(iniptb),2,quantile,0.95)
      tmpbias.ZouAhat = mean(tmpptb.ZouAhat)*tmpbias.ZouAhat%*%junkAhat2
      tmpbias.ZouAhat = c(tmpbias.ZouAhat*(1-tmpp0[not0]))
      tmptmpbias.ZouAhat = rep(0,pz)
      tmptmpbias.ZouAhat[not0] = tmpbias.ZouAhat
      tmpbias.ZouAhat = tmptmpbias.ZouAhat
    }		
    junkbias2.ZouAhat= rbind(junkbias2.ZouAhat,tmpbias.ZouAhat)  
    
    #### Uses tmpsd2 ######   
    tmpCI = cbind(tmpbhat-1.96*tmpsd2, tmpbhat+1.96*tmpsd2)                           #NormCovp, NormCovp Mean, NormLen
    tmpcov= (tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true)
    tmpcov.mean= (tmpCI[,1]<=beta.mean)*(tmpCI[,2]>=beta.mean)
    junklen.norm = rbind(junklen.norm, tmpCI[,2]-tmpCI[,1])
    junkcov.norm = rbind(junkcov.norm, tmpcov)
    junkcov.norm.mean = rbind(junkcov.norm.mean, tmpcov.mean)
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.norm=rbind(junkcov0.norm,tmpcov0)
    
    tmpCI = tmpCI + tmpbias.non0 #####added non0!
    tmpcov= (tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true)                                       #NormBCCovp, NormBCCovp Mean
    tmpcov.mean= (tmpCI[,1]<=beta.mean)*(tmpCI[,2]>=beta.mean)
    junkcov.normbc = rbind(junkcov.normbc, tmpcov) 
    junkcov.normbc.mean = rbind(junkcov.normbc.mean, tmpcov.mean) 
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.normbc=rbind(junkcov0.normbc,tmpcov0)
    
    #### Uses tmpsd ######   
    tmpCI = cbind(tmpbhat-1.96*tmpsd, tmpbhat+1.96*tmpsd)                           #NormCovp, NormCovp Mean, NormLen
    tmpcov= (tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true)
    tmpcov.mean= (tmpCI[,1]<=beta.mean)*(tmpCI[,2]>=beta.mean)
    junklen.norm2 = rbind(junklen.norm2, tmpCI[,2]-tmpCI[,1])
    junkcov.norm2 = rbind(junkcov.norm2, tmpcov)
    junkcov.norm.mean2 = rbind(junkcov.norm.mean2, tmpcov.mean)
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.norm2=rbind(junkcov0.norm2,tmpcov0) 
    
    tmpCI = tmpCI + tmpbias.non0
    tmpcov= (tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true)                                       #NormBCCovp, NormBCCovp Mean
    tmpcov.mean= (tmpCI[,1]<=beta.mean)*(tmpCI[,2]>=beta.mean)
    junkcov.normbc2 = rbind(junkcov.normbc2, tmpcov) 
    junkcov.normbc.mean2 = rbind(junkcov.normbc.mean2, tmpcov.mean) 
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.normbc2=rbind(junkcov0.normbc2,tmpcov0)
    
    tmpCI = cbind(tmpbhat-1.96*tmpsd, tmpbhat+1.96*tmpsd)
    tmpCI = tmpCI + tmpbias.nop0.non0
    tmpcov= (tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true)                                       #NormBCCovp, NormBCCovp Mean
    junkcov.normbc2.nop0 = rbind(junkcov.normbc2.nop0, tmpcov)  
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.normbc2.nop0=rbind(junkcov0.normbc2.nop0,tmpcov0)
    
    
    #### Uses tmpsd2.nop0 ######   
    tmpCI = cbind(tmpbhat-1.96*tmpsd2.nop0, tmpbhat+1.96*tmpsd2.nop0)                 #NormCovp, NormCovp Mean, NormLen
    tmpcov= (tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true)
    junklen.norm.nop0 = rbind(junklen.norm.nop0, tmpCI[,2]-tmpCI[,1])
    junkcov.norm.nop0 = rbind(junkcov.norm.nop0, tmpcov)
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.norm.nop0=rbind(junkcov0.norm.nop0,tmpcov0)
    
    tmpCI = tmpCI + tmpbias.nop0.non0 #####added non0!
    tmpcov= (tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true)                                       #NormBCCovp, NormBCCovp Mean
    junkcov.normbc.nop0 = rbind(junkcov.normbc, tmpcov) 
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.normbc.nop0=rbind(junkcov0.normbc.nop0,tmpcov0)
    
    #### Uses Zou2006 Ahat tmpsd.ZouAhat ######   
    tmpCI = cbind(tmpbhat-1.96*tmpsd.ZouAhat, tmpbhat+1.96*tmpsd.ZouAhat)                   #NormCovp, NormCovp Mean, NormLen
    tmpcov= (tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true)
    tmpcov.mean= (tmpCI[,1]<=beta.mean)*(tmpCI[,2]>=beta.mean)
    junklen.normA.ZouAhat = rbind(junklen.normA.ZouAhat, tmpCI[,2]-tmpCI[,1])
    junkcov.normA.ZouAhat = rbind(junkcov.normA.ZouAhat, tmpcov)
    junkcov.normA.mean.ZouAhat = rbind(junkcov.normA.mean.ZouAhat, tmpcov.mean)
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.normA.ZouAhat =rbind(junkcov0.normA.ZouAhat,tmpcov0)
    
    #### Quantiles of Betahatstar - betahat ####
    tmpCI = CI.quantlen(bptb=tmpptb[,-1],bhat=tmpbhat)
    junkAlp.opt = rbind(junkAlp.opt, as.vector(tmpCI[,1]))
    tmpcov = 1*((tmpCI[,2]<=beta.true)*(tmpCI[,3]>=beta.true))
    tmpcov.95 = 1*((tmpCI[,4]<=beta.true)*(tmpCI[,5]>=beta.true))
    tmpcov.mean = 1*((tmpCI[,2]<=beta.mean)*(tmpCI[,3]>=beta.mean))
    tmpcov.95.mean = 1*((tmpCI[,4]<=beta.mean)*(tmpCI[,5]>=beta.mean))
    tmpcov0 = 1-1*((tmpCI[,2]<=vec0)*(tmpCI[,3]>=vec0))
    tmpcov.950 = 1-1*((tmpCI[,4]<=vec0)*(tmpCI[,5]>=vec0))
    junkcov.quantlen = rbind(junkcov.quantlen, tmpcov)
    junkcov.95 = rbind(junkcov.95,tmpcov.95)
    junkcov.quantlen.mean = rbind(junkcov.quantlen.mean,tmpcov.mean)
    junkcov.95.mean = rbind(junkcov.95.mean, tmpcov.95.mean)
    junkcov.quantlen0 = rbind(junkcov.quantlen0, tmpcov0) 
    junkcov.950 = rbind(junkcov.950, tmpcov.950)
    junklen.95 = rbind(junklen.95, (tmpCI[,5] - tmpCI[,4]))
    junklen.quantlen = rbind(junklen.quantlen, (tmpCI[,3] - tmpCI[,2]))
    
    tmpCI = tmpCI + tmpbias
    tmpcov = 1*((tmpCI[,2]<=beta.true)*(tmpCI[,3]>=beta.true))
    tmpcov.95 = 1*((tmpCI[,4]<=beta.true)*(tmpCI[,5]>=beta.true))
    tmpcov0 = 1-1*((tmpCI[,2]<=vec0)*(tmpCI[,3]>=vec0))
    tmpcov.950 = 1-1*((tmpCI[,4]<=vec0)*(tmpCI[,5]>=vec0))
    junkcov.quantlenBC = rbind(junkcov.quantlen, tmpcov)
    junkcov.95BC = rbind(junkcov.95,tmpcov.95)
    junkcov.quantlen0BC = rbind(junkcov.quantlen0, tmpcov0) 
    junkcov.950BC = rbind(junkcov.950, tmpcov.950)
    
    #### Quantiles of betahatstar #####
    tmpCI = CI.quantlen(bptb=tmpptb[,-1])
    junkAlp.opt.star = rbind(junkAlp.opt.star, as.vector(tmpCI[,1]))
    tmpcov = 1*((tmpCI[,2]<=beta.true)*(tmpCI[,3]>=beta.true))
    tmpcov.95 = 1*((tmpCI[,4]<=beta.true)*(tmpCI[,5]>=beta.true))
    tmpcov.mean = 1*((tmpCI[,2]<=beta.mean)*(tmpCI[,3]>=beta.mean))
    tmpcov.95.mean = 1*((tmpCI[,4]<=beta.mean)*(tmpCI[,5]>=beta.mean))
    tmpcov0 = 1-1*((tmpCI[,2]<=vec0)*(tmpCI[,3]>=vec0))
    tmpcov.950 = 1-1*((tmpCI[,4]<=vec0)*(tmpCI[,5]>=vec0))
    junkcov.quantlen.star = rbind(junkcov.quantlen.star, tmpcov)
    junkcov.95.star = rbind(junkcov.95.star,tmpcov.95)
    junkcov.quantlen.mean.star = rbind(junkcov.quantlen.mean.star,tmpcov.mean)
    junkcov.95.mean.star = rbind(junkcov.95.mean.star, tmpcov.95.mean)
    junkcov.quantlen0.star = rbind(junkcov.quantlen0.star, tmpcov0) 
    junkcov.950.star = rbind(junkcov.950.star, tmpcov.950)
    junklen.95.star = rbind(junklen.95.star, (tmpCI[,5] - tmpCI[,4]))
    junklen.quantlen.star = rbind(junklen.quantlen.star, (tmpCI[,3] - tmpCI[,2]))
    
    aa=1
    
    p0.lo0 = pmin((sqrt(2/pi)*exp(-(aa*junklam[kk])/(4*mysig0^2))),0.49)
    p0.hi0 = pmax(0.5,pmin((1-(sqrt(2/pi)*exp(-(aa*junklam[kk])/(mysig0^2)))),.95))
    p0.lo = pmin((sqrt(2/pi)*exp(-(aa*junklam[kk])/(4*tmpsighat2))),0.49)
    p0.hi = pmax(0.5,pmin((1-(sqrt(2/pi)*exp(-(aa*junklam[kk])/(tmpsighat2)))),.95))
    print(c(p0.lo, p0.lo0, p0.hi, p0.hi0))
    
    tmpp0hi = c(1-1/sampsize^(0.25*(1+pz/sampsize)),pmin(1-1/sampsize^(0.25*(1+pz/sampsize)),0.95),
                (1-(sqrt(2/pi)*exp(-(aa*junklam[kk])/(4*mysig0^2)))), pmin((1-(sqrt(2/pi)*exp(-(aa*junklam[kk])/(4*mysig0^2)))),.95),
                (1-(sqrt(2/pi)*exp(-(aa*junklam[kk])/(mysig0^2)))),pmin((1-(sqrt(2/pi)*exp(-(aa*junklam[kk])/(mysig0^2)))),.95))
    tmpp0hi = c(tmpp0hi,(1-(sqrt(2/pi)*exp(-(aa*junklam[kk])/(4*tmpsighat2)))), pmin((1-(sqrt(2/pi)*exp(-(aa*junklam[kk])/(4*tmpsighat2)))),.95),
                (1-(sqrt(2/pi)*exp(-(aa*junklam[kk])/(tmpsighat2)))),pmin((1-(sqrt(2/pi)*exp(-(aa*junklam[kk])/(tmpsighat2)))),.95))
    junk.p0hi = rbind(junk.p0hi,tmpp0hi)
    
    tmpp0lo = c((log(sampsize)/sampsize)^(0.25*(1-pz/sampsize)),pmin((log(sampsize)/sampsize)^(0.25*(1-pz/sampsize)), 0.49),
                (sqrt(2/pi)*exp(-(aa*junklam[kk])/(mysig0^2))),pmin((sqrt(2/pi)*exp(-(aa*junklam[kk])/(mysig0^2))),0.49),
                (sqrt(2/pi)*exp(-(aa*junklam[kk])/(4*mysig0^2))), pmin((sqrt(2/pi)*exp(-(aa*junklam[kk])/(4*mysig0^2))),0.49))
    tmpp0lo = c(tmpp0lo,(sqrt(2/pi)*exp(-(aa*junklam[kk])/(tmpsighat2))),pmin((sqrt(2/pi)*exp(-(aa*junklam[kk])/(tmpsighat2))),0.49),
                (sqrt(2/pi)*exp(-(aa*junklam[kk])/(4*tmpsighat2))), pmin((sqrt(2/pi)*exp(-(aa*junklam[kk])/(4*tmpsighat2))),0.49))
    junk.p0lo = rbind(junk.p0lo,tmpp0lo)
    
    tmpp0.out = rbind("p0hi"=tmpp0hi,"p0lo"=tmpp0lo)
    colnames(tmpp0.out) = c("old","old.min","bound0","bound.min0","used0","used.min0","bound","bound.min","used","used.min")
    print(round(tmpp0.out,5))
    
    
    ##### HDR confidence intervals THE OLD WAY 
    tmpCImat.hdr0 = apply(rbind(tmpbias,beta.true, tmpptb[,-1]), 2, myhdrfun.ci.new,biasshift=FALSE)
    junkcov.hdr0 = rbind(junkcov.hdr0, tmpCImat.hdr0[1,])
    junkcov0.hdr0 = rbind(junkcov0.hdr0, tmpCImat.hdr0[2,])
    junklen.hdr0 = rbind(junklen.hdr0, tmpCImat.hdr0[3,])
    failurecount.ci0 = sum(tmpCImat.hdr0[4,]) + failurecount.ci0
    failurecount.ci.NA0 = sum(tmpCImat.hdr0[5,]) + failurecount.ci.NA0
    
    ##### HDR shifted
    tmpCImat.hdrBC = apply(rbind(tmpbias,beta.true,tmpptb[,-1]), 2, myhdrfun.ci.new,biasshift=TRUE)
    junkcov.hdrBC = rbind(junkcov.hdrBC, tmpCImat.hdrBC[1,])
    junkcov0.hdrBC = rbind(junkcov0.hdrBC, tmpCImat.hdrBC[2,])
    junklen.hdrBC = rbind(junklen.hdrBC, tmpCImat.hdrBC[3,])
    
    ##### SIMULTANEOUS CONFIDENCE INTERVALS
    tmpCImat.simul = simulci(rbind(beta.true,tmpbhat,tmpsd,tmpptb[,-1]))
    junkcov.simul = rbind(junkcov.simul,tmpCImat.simul[1])
    print(c(kk, tmpCImat.simul))
    junklen.simul = rbind(junklen.simul,tmpCImat.simul[2])
    junkgam.simul = rbind(junkgam.simul,tmpCImat.simul[3])
    
    print(c(kk,tmpCImat.simul))
    
    
    tmpCImat.simul = simulci(rbind(beta.true,tmpbhat.ini,tmpsd.ols,iniptb),type="OLS")
    junkcov.simul.ols = rbind(junkcov.simul.ols,tmpCImat.simul[1])
    junklen.simul.ols = rbind(junklen.simul.ols,tmpCImat.simul[2])
    junkgam.simul.ols = rbind(junkgam.simul.ols,tmpCImat.simul[3])
    
    
    ## USING BINI as a center
    #btrue,bhat,bhat.ini,bsig,bsig.ini,bptb,bptb.ini
    tmpCImat.simul = simulci.olsgamma(beta.true,tmpbhat,tmpbhat.ini,tmpsd,tmpsd.ols,tmpptb[,-1],iniptb)
    junkcov.simul2 = rbind(junkcov.simul2,tmpCImat.simul[1])
    junklen.simul2 = rbind(junklen.simul2,tmpCImat.simul[2])
    junkgam.simul2 = rbind(junkgam.simul2,tmpCImat.simul[3])
    
    #STABILIZE WITH SIGHAT/2
    #simulci <- function(bptb,type="Resample",sighat=0,nn=100) {
    tmpCImat.simul = simulci(rbind(beta.true,tmpbhat,tmpsd,tmpptb[,-1]),sighat=sqrt(tmpsighat2),nn=sampsize)
    junkcov.simul3 = rbind(junkcov.simul3,tmpCImat.simul[1])
    junklen.simul3 = rbind(junklen.simul3,tmpCImat.simul[2])
    junkgam.simul3 = rbind(junkgam.simul3,tmpCImat.simul[3])
    
    #SIMULTANEOUS HDR
    simul.alp = simulci.alp(tmpbhat,tmpsd,tmpptb[,-1])
    junksimul.alp = rbind(junksimul.alp,simul.alp)
    tmpCImat.hdr.simul = apply(rbind(tmpbias,beta.true, tmpptb[,-1]), 2, myhdrfun.ci.new,biasshift=FALSE,alp=simul.alp)
    junkcov.hdr.simul = rbind(junkcov.hdr.simul, prod(tmpCImat.hdr.simul[1,]))
    junklen.hdr.simul = rbind(junklen.hdr.simul, sum(tmpCImat.hdr.simul[3,]))
    
    #### Normal Uses tmpsd, c(0,0) if higher than p0.hi ######   
    tmpCI = cbind(tmpbhat-1.96*tmpsd, tmpbhat+1.96*tmpsd)                           #NormCovp, NormCovp Mean, NormLen
    tmpCI[tmpp0>p0.hi,] = c(0,0); #print(1*(tmpp0>p0.hi))
    tmpcov= (tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true)
    junklen.norm2p0 = rbind(junklen.norm2p0, tmpCI[,2]-tmpCI[,1])
    junkcov.norm2p0 = rbind(junkcov.norm2p0, tmpcov)
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.norm2p0=rbind(junkcov0.norm2p0,tmpcov0)
    
    #### Normal to parallel HDR when tp0 < plo
    tmpCI = t(apply(rbind(tmpbhat,tmpsd,tmpptb[,-1]), 2,mynormfun.ci))
    tmpcov= (tmpCI[,1]<=beta.true)*(tmpCI[,2]>=beta.true)
    junklen.norm2infl = rbind(junklen.norm2infl, tmpCI[,2]-tmpCI[,1])
    junkcov.norm2infl = rbind(junkcov.norm2infl, tmpcov)
    tmpcov0= 1 - 1*(tmpCI[,1]<=vec0)*(tmpCI[,2]>=vec0)
    junkcov0.norm2infl=rbind(junkcov0.norm2infl,tmpcov0)
    
    if (kk%%100==0) print(round(c(kk,apply(junkcov.hdr0,2,mean)),3))
  } 
  
  if(pz==20) {
    junk.p0lo = cbind(junk.p0lo, junk.p0lo)
    junk.p0hi = cbind(junk.p0hi, junk.p0hi)
  }
  result = rbind(
    "TrueBeta"				=	beta.true,                                          #TrueBeta 
    "EmpMean" 				= 	apply(junkhat,2,mean),                              #EmpMean
    "EmpMeanBC P/B^2"		=	apply(junkhat+junkbias2.nop0.non0,2,mean),          #EmpMeanBC P/B^2 with no P0 terms
    "EmpMeanBC ZouAhat"		=	apply(junkhat+junkbias2.ZouAhat,2,mean),            #EmpMeanBC ZouAhat
    "TypeIErr Empp0"		=	apply(1*(junkp0<0.05),2,mean),
    "ModelSize"				=	mean(junkmodel),                               #AveModelSize
    "EmpOLSSD"				=	apply(junkini,2,sd),                                       #EmpOLSSD 
    "AveOLSSD"				=	apply(junksd.ols,2,mean),                                  #AveOLSSD 
    "EmpSD"					=	apply(junkhat,2,sd),                                       #EmpSD 
    "AveSD (tmpsd)"			=	apply(junksd,2,mean),                                      #AveSD (tmpsd)
    "EmpBCSD"				=	apply(junkhat+junkbias2.non0,2,sd),                             #EmpBCSD
    "AveBCSD (tmpsd2)"		=	apply(junksd2,2,mean),                                     #AveBCSD (tmpsd2) 
    "AveBCSD (tmpsd2.nop0)" =	apply(junksd2.nop0,2,mean),
    "ZouAhatSD"				=	apply(junksd.ZouAhat,2,mean),                               #Zou2006 Ahat
    "EmpP0"					=	apply(junkhat==0,2,mean),                                  #EmpP0
    "AveP0"					=	apply(junkp0,2,mean),                                      #AveP0     
    "OLSCovp"				=	apply(junkcov.ols,2,mean),                                 #OLSCovp     
    "NormCovp tmpsd"		=	apply(junkcov.norm2,2,mean),                               #NormCovp tmpsd
    "NormCovp tmpsd p0hi"	=	apply(junkcov.norm2p0,2,mean),                               #NormCovp tmpsd
    "NormCovp tmpsd infl"	=	apply(junkcov.norm2infl,2,mean),
    "NormCovp tmpsd2"		=	apply(junkcov.norm,2,mean),                     #NormCovp tmpsd2
    "NormCovp tmpsd2.noP0"	=	apply(junkcov.norm.nop0,2,mean),                     #NormCovp tmpsd2
    "NormCovp ZouAhat" 		=	apply(junkcov.normA.ZouAhat,2,mean),                       #NormCovp ZouAhat
    "NormBCCovp tmpsd"		=	apply(junkcov.normbc2,2,mean),                             #NormBCCovp tmpsd
    "NormBCCovp tmpsd2"		=	apply(junkcov.normbc,2,mean),                              #NormBCCovp tmpsd2 
    "NormBCCovp(nop0) tmpsd2.nop0"		=    apply(junkcov.normbc.nop0,2,mean),        #NormBCCovp tmpsd2
    "NormBCCovp(nop0) tmpsd"		=	apply(junkcov.normbc2.nop0,2,mean),                             #NormBCCovp tmpsd 
    "HDRCovp"			=	apply(junkcov.hdr0,2,mean),
    "HDRCovpBC"				=	apply(junkcov.hdrBC,2,mean),
    "95QuantileCovp"		=	apply(junkcov.95,2,mean),                                  #95QuantileCovp
    "95QuantileCovp Star"	=	apply(junkcov.95.star,2,mean),                             #95QuantileCovp Star
    "SimulCovp"				=	mean(junkcov.simul),
    "SimulCovp GamOLS"		=	mean(junkcov.simul2),
    "SimulCovp GamSig"		=	mean(junkcov.simul3),
    "SimulCovp OLS"			=	mean(junkcov.simul.ols),
    "SimulCovp HDR"			=	mean(junkcov.hdr.simul),
    "OLSCovp0"				=	apply(junkcov0.ols,2,mean),                                #OLSCovp0
    "NormCovp0 tmpsd"		=	apply(junkcov0.norm2,2,mean),                              #NormCovp0 tmpsd
    "NormCovp0 tmpsd p0hi"	=	apply(junkcov0.norm2p0,2,mean),                              #NormCovp0 tmpsd
    "NormCovp0 tmpsd infl"	=	apply(junkcov0.norm2infl,2,mean),        
    "NormCovp0 tmpsd2"		=	apply(junkcov0.norm,2,mean),                              #NormCovp0 tmpsd2
    "NormCovp0 tmpsd2.nop0"		=	apply(junkcov0.norm.nop0,2,mean),                            #NormCovp0 tmpsd2
    "NormCovp0 ZouAhat"		=	apply(junkcov0.normA.ZouAhat,2,mean),                      #NormCovp0 ZouAhat
    "NormBCCovp0 tmpsd"		=	apply(junkcov0.normbc2,2,mean),                            #NormBCCovp0 tmpsd
    "NormBCCovp0 tmpsd2"	=	apply(junkcov0.normbc,2,mean),                             #NormBCCovp0 tmpsd2
    "NormBCCovp0(nop0) tmpsd2.nop0"		=    apply(junkcov0.normbc.nop0,2,mean),        #NormBCCovp tmpsd2
    "NormBCCovp0(nop0) tmpsd"		=	apply(junkcov0.normbc2.nop0,2,mean),                             #NormBCCovp tmpsd 
    "HDRCovp0"			=	apply(junkcov0.hdr0,2,mean),                                 #HDRCovp0
    "HDRCovp0BC"				=	apply(junkcov0.hdrBC,2,mean),
    "95QuantileCovp0" 		=	apply(junkcov.950,2,mean),                            #95QuantileCovp0
    "95QuantileCovp0 Star" 	=	apply(junkcov.950.star,2,mean),                   #95QuantileCovp0 Star
    "OLSLen"				=	apply(junklen.ols,2,mean),                                 #OLSLen
    "NormLen tmpsd"			=	apply(junklen.norm2,2,mean),                               #NormLen tmpsd
    "NormLen tmpsd p0hi"	=	apply(junklen.norm2p0,2,mean),                               #NormLen tmpsd
    "NormLen tmpsd infl"	=	apply(junklen.norm2infl,2,mean),
    "NormLen tmpsd2"		=	apply(junklen.norm,2,mean),                                #NormLen tmpsd2
    "NormLen tmpsd2.nop0"		=	apply(junklen.norm.nop0,2,mean),                                #NormLen tmpsd2
    "NormLen ZouAhat"		=	apply(junklen.normA.ZouAhat,2,mean),                       #NormLen ZouAhat
    "95QuantileLen"			=	apply(junklen.95,2,mean),
    "95QuantileLen Star"	=	apply(junklen.95.star,2,mean),                              #95QuantileLen Star
    "HDRLen"			=	apply(junklen.hdr0,2,mean),
    "HDRLenBC"				=	apply(junklen.hdrBC,2,mean),
    "SimulLen"				= mean(junklen.simul),
    "SimulLen GamOLS"		= mean(junklen.simul2),
    "SimulLen GamSig"		= mean(junklen.simul3),
    "SimulLen OLS"			= mean(junklen.simul.ols),
    "SimulLen HDR"			= mean(junklen.hdr.simul),
    "SimulGamma"			= mean(junkgam.simul),
    "SimulGamma GamOLS"		= mean(junkgam.simul2),
    "SimulGamma GamSig"		= mean(junkgam.simul3),
    "SimulGamma OLS"		= mean(junkgam.simul.ols),
    "SimulAlpha"			= mean(junksimul.alp),
    "Lambda"				= mean(junklam),
    "LambdaPtb"				=	mean(junklam.ptb),
    "Sighat2"				=	mean(junksig2),
    "junk.p0lo"				= apply(junk.p0lo,2,mean),
    "junk.p0hi"				= apply(junk.p0hi,2,mean)
  )
  
  
  
  colnames(result)=paste("X",1:pz,sep="")
  print(round(result,4))
  
  (nptb0 <- length(junk[,1]))
  
  
  
  if(cluster == FALSE) setwd(outputwd)
  write.csv(round(result,4),paste(outputwd.sum,"summary.b",sprintf("%1d",beta.true[4]*10),"n",sampsize,"ptb",nptb0,"pz",pz,"sig",mysig0,"corX",cor.X,".BIC.csv",sep=""))
  
  par(mfrow=c(1,2))
  plot(density(junkgam.simul))
  plot(density(junkgam.simul.ols))
  par(mfrow=c(1,2))
  hist(junkgam.simul)
  hist(junkgam.simul.ols)
  #dev.off()
  if(cluster == FALSE) setwd(inputwd)
  
  
  paste("summary.",sampsize,"ptb",nptb0,"pz",pz,"b",sprintf("%1d",beta.true[4]*10),"sig",mysig0,"corX",cor.X,".csv",sep="")
  
  
  print(failurecount.ci)
  print(failurecount.ci.NA)
  print(failurecount.mode)
  print(failurecount.ci0)
  print(failurecount.ci.NA0)
  print(failurecount.mode0)
  
  print(round(apply(junk.p0lo,2,mean),5))
  print(round(apply(junk.p0hi,2,mean),5))
  
  
}else{print(paste(tmpcsv, "exists"))}                                                              
