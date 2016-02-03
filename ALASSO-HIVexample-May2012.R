# Author: Jessica Minnier
# This code runs the analysis of the Stanford HIV Data described in the article. This can be adapted to run similar analyses on other data sets.


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


myseed = 103


inputwd = getwd()
outputwd = inputwd
sourcewd = inputwd
outputwdtex = outputwd


setwd(sourcewd)
source("functions-regression.R")
source("functions-CIs.R")


################################################################################################################
############## Ahat is Ahat + diag P(betaptb!=0)/ols est (from Tianxi's paper)
################################################################################################################

setwd(inputwd)
hivdata = read.table("PI_DATA.cleaned.txt",header=TRUE)
mydata0 = as.matrix(hivdata)
mydata0[,1] = log(mydata0[,1])

feat0 = mydata0[,-1]
ind.rm = apply(feat0, 2, mean, na.rm=T)
feat.all0 = feat0[,ind.rm>=.005] #removes mutations with less than .5% prevalence, 15 mutations
mu.cov = apply(feat.all0,2,mean,na.rm=T)
#sd.cov = apply(mydata0[,-1],2,sd,na.rm=T)
#feat.st = (mydata0[,-1] - VTM(mu.cov, nrow(mydata0)))/VTM(sd.cov,nrow(mydata0))
feat.st = feat.all0 -VTM(mu.cov, nrow(mydata0))
mydata = cbind(mydata0[,1] - mean(mydata0[,1]),feat.st)
mydata.unst = cbind(mydata0[,1],feat.all0)

mynames = names(hivdata)[-1]
mynames = mynames[ind.rm>=0.005]
mynames = c(names(hivdata)[1],mynames)
#############################
nptb = 500
tmpn = sampsize = dim(mydata)[1]
pz = dim(mydata)[2] - 1
#############################

set.seed(myseed)

beta.hat = Est.ALASSO(mydata,bictype="min")
beta.ptb.all = NULL
for(i in 1:nptb)
{
  if(i%%100==0) print(i)
  tmpG = rexp(nrow(mydata)) 
  means2 = apply(mydata.unst*tmpG,2,mean)
  mydata.c2 = (mydata.unst - VTM(means2,sampsize))*sqrt(tmpG)
  beta.ptb = Est.ALASSO(mydata.c2,rtn="notall",s0=beta.hat[[1]][1],bictype="min") #don't return ptb'd Ahat
  beta.ptb.all = rbind(beta.ptb.all,beta.ptb)
}

output = list("beta.out" = beta.hat$beta.out,"Ahat" = beta.hat$Ahat, "beta.ptb.all" = beta.ptb.all)


###########################################################################################
#############   Summary    ############################################################
###########################################################################################


junk = c(beta.hat$beta.out, beta.hat$Ahat, beta.ptb.all)
junk = as.matrix(junk);mode(junk)="numeric"

junklam  = beta.hat$beta.out[1]; 
junkhat  = beta.hat$beta.out[1+1:pz]; 			#alasso estimates
junkini  = beta.hat$beta.out[1+pz+1:pz];		#ols estimates
junkbias = beta.hat$beta.out[1+2*pz+1:pz]; 		#bias estimate
junkAhat = c(beta.hat$Ahat)		#X'X/n
junksighat2 = beta.hat$sighat2
junkptb  = beta.ptb.all	#ptbed

beta.mean <-mean(junkhat)
vec0 <- rep(0,pz)

############# in summary, the loop through data sets starts here #############


tmpptb = junkptb
iniptb = tmpptb[,-(1:(1+pz))]; tmpptb = tmpptb[,1:(1+pz)]

tmpbhat = junkhat
tmpmodel = length(tmpbhat[tmpbhat!=0])		                                                #model size
tmpbhat.ini = as.numeric(junkini)
tmpsighat2 = junksighat2


dvec = apply(tmpptb!=0,2,mean)[-1]/junkini^2                                           #Ahat = Ahat + p(b*!=0)/bhat^2
junkAhat = junkAhat + (junklam/tmpn)*as.vector(diag(dvec))


not0 = which(tmpbhat!=0)				                                                    #Zou2006 Ahat (not Ahat but Jacobian)
Ahat.submat = tmpn*matrix(junkAhat,nrow=pz)[not0,not0]
mat1 = solve(Ahat.submat + (junklam/tmpn)*diag(1/(abs(tmpbhat*junkini)[not0])))
#######need estimate of mysig0
fit = lm(mydata[,1]~mydata[,-1])
mysig0 = (predict(fit,se.fit=T))$res

junkAhat2 = mysig0*(mat1%*%Ahat.submat%*%mat1)
tmpsd.ZouAhat = rep(0,pz)
tmpsd.ZouAhat[not0] = sqrt(diag(junkAhat2))

tmpsd.naive = naivesd(tmpbhat, junkAhat,tmpn)


############ emp sd of ptb betas #############################################################

tmpsd = sqrt(apply( (tmpptb[,-1]-VTM(junkhat,nrow(tmpptb)))^2,2,mean));                #tmpsd: 

tmpsd.ols=apply(iniptb,2,sd)                                             #OLSCovp, OLSCovp Mean, OLSLen
tmpCI.ols=cbind(junkini-1.96*apply(iniptb,2,sd),junkini+1.96*apply(iniptb,2,sd))      #empirical sd of iniptb
tmplen.ols=tmpCI.ols[,2]-tmpCI.ols[,1]
tmpcov0.ols = 1 - 1*(tmpCI.ols[,1]<=vec0)*(tmpCI.ols[,2]>=vec0)


tmpp0 = apply(tmpptb[,-1]==0,2,mean)			                         #P(beta=0)

tmpbias=(-1)^(apply(sign(tmpptb[,-1])==1,2,mean)<apply(sign(tmpptb[,-1])==-1,2,mean))/
  max(abs(quantile(iniptb,c(0.075,0.975))))                                     #apply(abs(iniptb),2,quantile,0.95)
tmpbias = mean(tmpptb[,1])*tmpbias%*%solve(matrix(junkAhat,ncol=pz))/tmpn
tmpbias = c(tmpbias*(1-tmpp0)) 				       
tmpbias.non0 = tmpbias*(tmpbhat!=0)                                         #changed from 1-tmpp0

tmpsd2 = sqrt(tmpsd^2+tmpbias^2)

### EmpMeanBC P/B
tmpbhatbc = tmpbhat + tmpbias.non0

#### Uses tmpsd ######   
#p0.hi = pmin(1-1/sampsize^(0.25*(1+pz/sampsize)),0.95)
aa=1
p0.lo = pmin((sqrt(2/pi)*exp(-(aa*junklam)/(4*tmpsighat2))),0.49)
p0.hi = pmax(0.5,pmin((1-(sqrt(2/pi)*exp(-(aa*junklam)/(tmpsighat2)))),.95))


tmpCI.norm = cbind(junkhat-1.96*tmpsd, junkhat+1.96*tmpsd)                           #NormCovp, NormCovp Mean, NormLen
tmpCI.norm[tmpp0>p0.hi,] = c(0,0)
tmplen.norm = tmpCI.norm[,2] - tmpCI.norm[,1]
tmpcov0.norm= 1 - 1*(tmpCI.norm[,1]<=vec0)*(tmpCI.norm[,2]>=vec0)

tmpCI.normbc = tmpCI.norm + tmpbias.non0
tmpcov0.normbc = 1 - 1*(tmpCI.normbc[,1]<=vec0)*(tmpCI.normbc[,2]>=vec0)


#### Uses tmpsd2 ######   

tmpCI.norm2 = cbind(junkhat-1.96*tmpsd2, junkhat+1.96*tmpsd2)                           #NormCovp, NormCovp Mean, NormLen
tmplen.norm2 = tmpCI.norm2[,2] - tmpCI.norm2[,1]
tmpcov0.norm2= 1 - 1*(tmpCI.norm2[,1]<=vec0)*(tmpCI.norm2[,2]>=vec0)

tmpCI.normbc2 = tmpCI.norm2 + tmpbias.non0
tmpcov0.normbc2 = 1 - 1*(tmpCI.normbc2[,1]<=vec0)*(tmpCI.normbc2[,2]>=vec0)


#### Quantiles of betahatstar #####

tmpCI.95quant = CI.quantlen(bptb=tmpptb[,-1])[,4:5]
tmpcov0.95quant = 1-1*((tmpCI.95quant[,1]<=vec0)*(tmpCI.95quant[,2]>=vec0))
tmplen.95quant = tmpCI.95quant[,2] - tmpCI.95quant[,1]


#### Uses Zou2006 Ahat tmpsd.ZouAhat ######   

#need estimate of mysig0
tmpCI.Zou = cbind(junkhat-1.96*tmpsd.ZouAhat, junkhat+1.96*tmpsd.ZouAhat)                   #NormCovp, NormCovp Mean, NormLen
tmplen.Zou = tmpCI.Zou[,2]-tmpCI.Zou[,1]
tmpcov0.Zou = 1 - 1*(tmpCI.Zou[,1]<=vec0)*(tmpCI.Zou[,2]>=vec0)


#  p0.lo =pmin((log(sampsize)/sampsize)^(0.25*(1-pz/sampsize)), 0.49)
#  p0.hi = pmin(1-1/sampsize^(0.25*(1+pz/sampsize)),0.95)

##### HDR confidence intervals
#change function to output actual HDR
tmpCImat.hdr0 = apply(rbind(tmpbias, rep(0,pz), tmpptb[,-1]), 2, myhdrfun.ci.new2,biasshift=FALSE,rtn="CI")

#for uneven lengths of HDRs, outputs list; fill list with NA's to make same dimensions
maxlen = max(sapply(tmpCImat.hdr0, length,simplify=TRUE))
stretch = function(k) {
  newvec = rep(NA,maxlen)
  newvec[1:length(k)] = k
  newvec
}
tmpCImat.hdr = sapply(tmpCImat.hdr0, stretch, simplify=T)

#maxlen = 7
#tmpCImat.hdr = tmpCImat.hdr0
tmpCI.hdr = t(tmpCImat.hdr[5:maxlen,]) #first column is 1 if the mass at 0 is included
tmplen.hdr = tmpCImat.hdr[3,]
tmpcov0.hdr = tmpCImat.hdr[4,]
failurecount.ci = sum(tmpCImat.hdr[1,])
failurecount.ci.NA = sum(tmpCImat.hdr[2,])

#simultaneous ci's

##### SIMULTANEOUS CONFIDENCE INTERVALS
tmpCImat.simul = simulci(rbind(rep(0,pz),tmpbhat,tmpsd,tmpptb[,-1]),rtn="Reg")
tmpcov.simul = tmpCImat.simul[1]
tmplen.simul = tmpCImat.simul[2]
tmpgam.simul = tmpCImat.simul[3]

tmpCImat.simul.ols = simulci(rbind(rep(0,pz),junkini,apply(iniptb,2,sd),iniptb),type="OLS")
tmpcov.simul.ols = tmpCImat.simul.ols[1]
tmplen.simul.ols = tmpCImat.simul.ols[2]
tmpgam.simul.ols = tmpCImat.simul.ols[3]

tmpCI.simul = t(simulci(rbind(rep(0,pz),tmpbhat,tmpsd,tmpptb[,-1]),rtn="CI"))

simul.alp = simulci.alp(tmpbhat,tmpsd,tmpptb[,-1])
tmpCImat.hdr.simul0 = apply(rbind(tmpbias,rep(0,pz),tmpptb[,-1]), 2, myhdrfun.ci.new,biasshift=FALSE,alp=simul.alp,rtn="CI")
maxlen = max(sapply(tmpCImat.hdr.simul0, length,simplify=TRUE))
stretch = function(k) {
  newvec = rep(NA,maxlen)
  newvec[1:length(k)] = k
  newvec
}
tmpCImat.hdr.simul = sapply(tmpCImat.hdr.simul0, stretch, simplify=T)

tmpcov.hdr.simul = prod(tmpCImat.hdr.simul[1,])
tmplen.hdr.simul = sum(tmpCImat.hdr.simul[3,])
tmpCI.hdr.simul = t(tmpCImat.hdr.simul[5:maxlen,])



tp0 = apply(tmpptb, 2, function(k) sum(k==0)/length(k))


result.estimates = rbind(
  "OLS" = junkini,
  "AdaptiveLasso" = tmpbhat,
  "AdaptiveLassoBC" = tmpbhatbc,
  "Prob=0" = tp0[-1],
  "MAF" = apply(mydata,2,min)[-1]
)
colnames(result.estimates) = mynames[-1]
round(result.estimates,4)

result.modelsize = rbind(
  "Total Number" = pz,
  "AdaptiveLasso" = tmpmodel
)
result.modelsize


result.ci = cbind(
  "OLS" = tmpCI.ols,
  "Norm tmpsd" = tmpCI.norm,
  "NormBC tmpsd" = tmpCI.normbc,
  "Norm tmpsd2" = tmpCI.norm2,
  "NormBC tmpsd2" = tmpCI.normbc2,
  "95Quantile" = tmpCI.95quant,
  "Norm Zou" = tmpCI.Zou,
  "HDR" = tmpCI.hdr
)
colnames(result.ci) = c(
  paste(rep(c("OLS","Norm tmpsd", "NormBC tmpsd","Norm tmpsd2", "NormBC tmpsd2", "95Quantile","Norm Zou"),each=2), 1:2),
  paste("HDR", 1:(dim(tmpCI.hdr))[2])
)
rownames(result.ci) = mynames[-1]
result.ci

result.ci.len = rbind(
  "OLS" = tmplen.ols,
  "Norm tmpsd" = tmplen.norm,
  "NormBC tmpsd" = tmplen.norm,
  "Norm tmpsd2" = tmplen.norm2,
  "NormBC tmpsd2" = tmplen.norm2,
  "95Quantile" = tmplen.95quant,
  "HDR" = tmplen.hdr,
  "Norm Zou" = tmplen.Zou,
  "Prob=0" = tp0[-1]
)
colnames(result.ci.len) = mynames[-1]
result.ci.len
apply(result.ci.len,1,mean)

result.ci.cov0 = rbind(
  "OLS" = tmpcov0.ols,
  "Norm tmpsd" = tmpcov0.norm,
  "NormBC tmpsd" = tmpcov0.normbc,
  "Norm tmpsd2" = tmpcov0.norm2,
  "NormBC tmpsd2" = tmpcov0.normbc2,
  "95Quantile" = tmpcov0.95quant,
  "HDR" = tmpcov0.hdr,
  "Norm Zou" = tmpcov0.Zou
)
colnames(result.ci.cov0) = mynames[-1]
result.ci.cov0
apply(result.ci.cov0,1,sum)

result.sd = rbind(
  "OLS sd" = tmpsd.ols,
  "tmpsd" = tmpsd,
  "tmpsd2" = tmpsd2,
  "Zou sd" = tmpsd.ZouAhat
)
colnames(result.sd) = mynames[-1]
result.sd

####################################################################################
# 	plots and tables							
####################################################################################
setwd(outputwdtex)

####################################################################################
# 	non0 beta plot							
####################################################################################

length(tmpbhatbc)
#non0 = which(tmpbhatbc!=0)
non0 = which(result.ci.cov0[2,]!=0)
bhatnon0 = tmpbhatbc[non0]
names(bhatnon0) = NULL
namesnon0 = mynames[-1]
namesnon0 = namesnon0[non0]

#pdf("betabarplot.pdf",width=10.5,height=6)
postscript("betabarplot.eps",width=10.5,height=6,horizontal=T)
par(oma=c(8,6,0,0))
#mp = barplot(bhatnon0,ylab = expression(hat(beta)[BC]),cex.lab=.8)
mp = barplot(bhatnon0,ylab = expression(hat(beta)),cex.lab=.8)

starsig = function(vec) {
  vecf = as.factor(vec)
  if(nlevels(vecf)==1) {levels(vecf) = c("*")}
  else{levels(vecf) = c(" ", "*")}
  as.character(vecf)}
xstart = -2.4
ystart = -1.4
incr = .23

textvec = c("Codon",namesnon0)
text(c(xstart,mp),ystart,textvec,xpd = NA,cex=1.2)
textvec = c("Johnson et al", starsig(c(1,0,1,0,0,1,1,0,1,1,0,0,1,1,0,0,1,0)))
text(c(xstart,mp),ystart - incr,textvec,xpd = NA,cex=1.2)
textvec = c("Wu Adjusted", starsig(c(1,1,1,1,0,1,1,1,1,1,0,0,1,1,1,0,1,0)))
text(c(xstart,mp),ystart - 2*incr,textvec,xpd = NA,cex=1.2)

#typevec = c("OLS","Zou","Perturbation (Quantile)","Perturbation (HDR)","Perturbation (Normal)")
#typeind = c(1,8,6,7,2)

typevec = c("OLS","Asymptotic Based","Perturbation (Normal)","Perturbation (HDR)","Perturbation (Q)")
typeind = c(1,8,2,7,6)
for(jj in 1:(length(typeind))) {
  textvec = c(typevec[jj],starsig(result.ci.cov0[typeind[jj],non0]))
  text(c(xstart,mp),(ystart - incr*(jj+2)),textvec,xpd = NA,cex=1.2)
}

textvec = c(round(tp0[non0+1],2))
text(xstart,(ystart - incr*(jj+3)), expression(hat(p)[0]),xpd=NA,cex=1.2)
#text(c(xstart,mp),(ystart - incr*(jj+3)),textvec,xpd=NA,cex=1.2)
text(mp,(ystart - incr*(jj+3)),textvec,xpd=NA,cex=.8)	

dev.off()


####################################################################################
# 	confidence interval plot NORMAL				
####################################################################################

#onem.ind = which(apply(mydata[,-1],2,table)[2,]<3)
#myci = tmpCI.norm[-onem.ind,]
#myest = tmpbhatbc[-onem.ind]


myci = tmpCI.norm
myest = tmpbhat

non0 = which(tmpbhatbc!=0)
cov0.wu = rep(0,length(tmpbhatbc))
cov0.wu[non0] = c(1,1,1,1,1,1,1,1,1,1,1,1,1)
colorR= -1*tmpcov0.norm +2*cov0.wu + 1 #black if both 0, red if both 1, green if only Wu
color= colorR
color = rep(1, dim(myci)[1])

linetype = rep(1, dim(myci)[1])
#linetype= -1*tmpcov0.norm +5*cov0.wu + 1 #1 if both 0, 5 if both 1, 6 if only non0 estimate

pdf("ci-normptb.pdf",width=7,height=10)
#postscript("ci-normptb.eps",width=7,height=10)
axislab <- c()
plot(c(min(myci[,1]),max(myci[,2])),c(0,dim(myci)[1]),type="n",xlab=expression(beta),ylab="",yaxt='n',main="95% Perturbation (Normal) CIs")
abline(v=0,lty=2)
for(ii in 1:dim(myci)[1]) {
  #axislab <- cbind(axislab,k)
  lines(myci[ii,],c(ii,ii),col=color[ii],lty=linetype[ii])
  points(myest[ii],ii,col=color[ii],pch=20)
  if((myci[ii,2]-myci[ii,1])==0) {points(0,ii,col=color[ii],pch=8,cex=1)}
}

#axis(2,1:dim(myci)[1],mynames[-1])
axis(2,which(colorR>1),mynames[(which(colorR>1)+1)],cex.axis=.6,las=T)
#legend("topright",inset=.05,c("R rq() Output", "MCMC Ptb: Frequentist", "MCMC Lapl: Bayesian Posterior Inference","True Beta"),col=c(color,"green"),lty=c(1,1,1,-1),pch=17,cex=.7)
dev.off()


####################################################################################
# 	confidence interval plot HDR				
####################################################################################

#myci = tmpCI.hdr[-onem.ind,]
#myest = tmpbhatbc[-onem.ind]

myci = tmpCI.hdr
myest = tmpbhat

non0 = which(tmpbhatbc!=0)
cov0.wu = rep(0,length(tmpbhatbc))
cov0.wu[non0] = 1 #c(1,1,1,1,1,1,1,1,1,1,1,1,1)
colorR= -1*tmpcov0.hdr +2*cov0.wu + 1 #black if both 0, red if both 1, green if only Wu
#color= colorR
color = rep(1, dim(myci)[1])



#linetype= -1*tmpcov0.hdr +5*cov0.wu + 1 #1 if both 0, 5 if both 1, 6 if only non0 estimate
linetype = rep(1, dim(myci)[1])

pdf("ci-hdrptb.pdf",width=7,height=10)
axislab <- c()
plot(c(min(myci,na.rm=T),max(myci[,2:dim(myci)[2]],na.rm=T)),c(0,dim(myci)[1]),type="n",xlab=expression(beta),ylab="",yaxt='n',main="95% Perturbation (HDR) CRs")

abline(v=0,lty=2)

for(ii in 1:dim(myci)[1]) {
  #axislab <- cbind(axislab,k)
  if(myci[ii,1]==1) {points(0,ii,col=color[ii],pch=8,cex=1)}
  lines(myci[ii,-1],rep(ii,length(myci[ii,])-1),col=color[ii],lty=linetype[ii])
  points(myest[ii],ii,col=color[ii],pch=20,cex=1)
}

#axis(2,1:dim(myci)[1],mynames[-1])
axis(2,which(colorR>1),mynames[(which(colorR>1)+1)],cex.axis=.6,las=T)
#legend("topright",inset=.05,c("R rq() Output", "MCMC Ptb: Frequentist", "MCMC Lapl: Bayesian Posterior Inference","True Beta"),col=c(color,"green"),lty=c(1,1,1,-1),pch=17,cex=.7)
dev.off()



####################################################################################
# 	confidence interval plot HDR and NORMAL	with SIMULTANEOUS CONF REGIONS
####################################################################################
shade <- function(xl,xr, yl, yu, ...)
{
  nl <- length(xl)
  for(i in seq(1,(nl-1),by=2))
  {
    polygon(c(xl[i], xr[i], xr[i+1], xl[i+1]),
            c(yl[i], yl[i], yu[i+1], yu[i+1]), ...)
  }
}


#pdf("ci-both.pdf",width=14,height=10)
postscript(file="ci-both.eps",horizontal=T,width=14,height=10)
par(mfrow=c(1,2))

#pdf("ci-both.pdf",width=7,height=10)

myci = tmpCI.norm
myest = tmpbhat

#non0 = which(tmpbhatbc!=0)
non0 = which(result.ci.cov0[2,]!=0)
cov0.wu = rep(0,length(tmpbhatbc))
#cov0.wu[non0] = c(1,1,1,1,1,1,1,1,1,1,1,1,1)
colorR= -1*tmpcov0.norm +2*cov0.wu + 1 #black if both 0, red if both 1, green if only Wu
color= colorR
color = rep(1, dim(myci)[1])

linetype = rep(1, dim(myci)[1])
#linetype= -1*tmpcov0.norm +5*cov0.wu + 1 #1 if both 0, 5 if both 1, 6 if only non0 estimate


axislab <- c()
plot(c(min(tmpCI.simul[,1]),max(tmpCI.simul[,2])),c(0,2*dim(myci)[1]),type="n",xlab=expression(beta),ylab="",yaxt='n',main="95% Perturbation (Normal) CIs")
abline(v=0,lty=2)
shade(rep(tmpCI.simul[,1],each=2),rep(tmpCI.simul[,2],each=2),1:(2*(dim(myci)[1]))+0.5,1:(2*(dim(myci)[1]))+0.5,col="lightgrey",border=NA)


for(ii in 1:dim(myci)[1]) {
  #axislab <- cbind(axislab,k)
  lines(myci[ii,],c(2*ii,2*ii),col=color[ii],lty=linetype[ii],lwd=2)
  points(myest[ii],2*ii,col=color[ii],pch=20)
  if((myci[ii,2]-myci[ii,1])==0) {points(0,2*ii,col=color[ii],pch=8,cex=1)}
}



#axis(2,1:dim(myci)[1],mynames[-1])
axis(2,2*non0,mynames[non0+1],cex.axis=.6,las=T)
#legend("topright",inset=.05,c("R rq() Output", "MCMC Ptb: Frequentist", "MCMC Lapl: Bayesian Posterior Inference","True Beta"),col=c(color,"green"),lty=c(1,1,1,-1),pch=17,cex=.7)

####********************************************************************************************************###
#### HDR
####********************************************************************************************************###

myci = tmpCI.hdr
myest = tmpbhat

#non0 = which(tmpbhatbc!=0)
non0 = which(result.ci.cov0[7,]!=0)
cov0.wu = rep(0,length(tmpbhatbc))
#cov0.wu[non0] = c(1,1,1,1,1,1,1,1,1,1,1,1,1)
colorR= -1*tmpcov0.hdr +2*cov0.wu + 1 #black if both 0, red if both 1, green if only Wu
#color= colorR
color = rep(1, dim(myci)[1])

#linetype= -1*tmpcov0.hdr +5*cov0.wu + 1 #1 if both 0, 5 if both 1, 6 if only non0 estimate
linetype = rep(1, dim(myci)[1])

axislab <- c()
plot(c(min(tmpCI.simul[,1],na.rm=T),max(tmpCI.simul[,2],na.rm=T)),c(0,2*dim(myci)[1]),type="n",xlab=expression(beta),ylab="",yaxt='n',main="95% Perturbation (HDR) CRs")

abline(v=0,lty=2)
shade(rep(tmpCI.hdr.simul[,2],each=2),rep(tmpCI.hdr.simul[,3],each=2),1:(2*(dim(myci)[1]))+0.5,1:(2*(dim(myci)[1]))+0.5,col="lightgrey",border=NA)
shade(rep(tmpCI.hdr.simul[,4],each=2),rep(tmpCI.hdr.simul[,5],each=2),1:(2*(dim(myci)[1]))+0.5,1:(2*(dim(myci)[1]))+0.5,col="lightgrey",border=NA)

for(ii in 1:dim(myci)[1]) {
  #axislab <- cbind(axislab,k)
  if(myci[ii,1]==1) {points(0,2*ii,col=color[ii],pch=8,cex=1)}
  if(tmpCI.hdr.simul[ii,1]==1){points(0,2*ii,col="lightgrey",pch=8,cex=1.5)}
  yval = rep(2*ii,length(myci[ii,])-1)
  xval = myci[ii,-1]
  if(length(xval)>2) {xval = c(xval[1:2],NA,xval[3:4]); yval = c(yval[1:2],NA,yval[3:4])}
  lines(xval,yval,col=color[ii],lty=linetype[ii],lwd=2)
  points(myest[ii],2*ii,col=color[ii],pch=20,cex=1)
}

#axis(2,1:dim(myci)[1],mynames[-1])
axis(2,2*non0,mynames[non0+1],cex.axis=.6,las=T)
#legend("topright",inset=.05,c("R rq() Output", "MCMC Ptb: Frequentist", "MCMC Lapl: Bayesian Posterior Inference","True Beta"),col=c(color,"green"),lty=c(1,1,1,-1),pch=17,cex=.7)


dev.off()

####################################################################################
#	COLOR VERSION: May 5, 2012
# 	confidence interval plot HDR and NORMAL	with SIMULTANEOUS CONF REGIONS
####################################################################################
shade <- function(xl,xr, yl, yu, ...)
{
  nl <- length(xl)
  for(i in seq(1,(nl-1),by=2))
  {
    polygon(c(xl[i], xr[i], xr[i+1], xl[i+1]),
            c(yl[i], yl[i], yu[i+1], yu[i+1]), ...)
  }
}

wuname = c(10,30,32,33,46,47,48,50,54,76,84,88,90) #exactly same as nonzer0
wuname = paste("P",wuname,sep="")

#pdf("ci-both.pdf",width=14,height=10)
#postscript(file="ci-both-color.eps",horizontal=T,width=14,height=10)
png(file="ci-both-color.png",width=1.4*480*1.5,height=480*1.5)
par(mfrow=c(1,2))

#NORMAL CR

myci0 = tmpCI.norm
myest0 = tmpbhat
mycisim0 = tmpCI.simul
cilen = -1*abs(myci0[,1]-myci0[,2])

tmpind = order(1*(myest0!=0),myest0,cilen)
myest = myest0[tmpind]
myci = myci0[tmpind,]
mycisim = mycisim0[tmpind,]
mynames2 = mynames[tmpind]

#non0 = which(tmpbhatbc!=0)
non0 = which(result.ci.cov0[2,tmpind]!=0)
cov0.wu = rep(0,length(tmpbhatbc))
cov0.wu[match(wuname,colnames(result.ci.cov0))] = 1 
colorR= 2*tmpcov0.norm[tmpind] -1*cov0.wu[tmpind] + 1 #black if both 0, red if both 1, green if only Wu
color= colorR
#color = rep(1, dim(myci)[1])

linetype = rep(1, dim(myci)[1])
#linetype= -1*tmpcov0.norm +5*cov0.wu + 1 #1 if both 0, 5 if both 1, 6 if only non0 estimate


axislab <- c()
plot(c(min(mycisim[,1]),max(mycisim[,2])),c(0,2*dim(myci)[1]),type="n",xlab=expression(beta),ylab="",yaxt='n',main="95% Perturbation (Normal) CIs")
abline(v=0,lty=2)
shade(rep(mycisim[,1],each=2),rep(mycisim[,2],each=2),1:(2*(dim(myci)[1]))+0.5,1:(2*(dim(myci)[1]))+0.5,col="lightgrey",border=NA)


for(ii in 1:dim(myci)[1]) {
  #axislab <- cbind(axislab,k)
  lines(myci[ii,],c(2*ii,2*ii),col=color[ii],lty=linetype[ii],lwd=2)
  points(myest[ii],2*ii,col=color[ii],pch=20)
  if((myci[ii,2]-myci[ii,1])==0) {points(0,2*ii,col=color[ii],pch=8,cex=1)}
}



#axis(2,1:dim(myci)[1],mynames[-1])
axis(2,2*non0,names(non0),cex.axis=.6,las=T)
#legend("topright",inset=.05,c("R rq() Output", "MCMC Ptb: Frequentist", "MCMC Lapl: Bayesian Posterior Inference","True Beta"),col=c(color,"green"),lty=c(1,1,1,-1),pch=17,cex=.7)



#### HDR

myci0 = tmpCI.hdr
myest0 = tmpbhat
mycisim0 = tmpCI.hdr.simul
cilen = (abs(myci0[,2]-myci0[,3]))
cilen2 = abs(myci0[,5]-myci0[,4])
cilen2[is.na(cilen2)] <- 0
cilen = cilen + cilen2
cilen = -1*cilen

tmpind = order(1*(myest0!=0),myest0,cilen)
myest = myest0[tmpind]
myci = myci0[tmpind,]
mycisim = mycisim0[tmpind,]
mynames2 = mynames[tmpind]

non0 = which(result.ci.cov0[7, tmpind]!=0)
cov0.wu = rep(0,length(tmpbhatbc))
cov0.wu[match(wuname,colnames(result.ci.cov0))] = 1 
colorR= 2*tmpcov0.norm[tmpind] -1*cov0.wu[tmpind] + 1 #black if both 0, red if both 1, green if only Wu
color= colorR
#color = rep(1, dim(myci)[1])

linetype = rep(1, dim(myci)[1])
#linetype= -1*tmpcov0.norm +5*cov0.wu + 1 #1 if both 0, 5 if both 1, 6 if only non0 estimate



axislab <- c()
plot(c(min(tmpCI.simul[,1],na.rm=T),max(tmpCI.simul[,2],na.rm=T)),c(0,2*dim(myci)[1]),type="n",xlab=expression(beta),ylab="",yaxt='n',main="95% Perturbation (HDR) CRs")

abline(v=0,lty=2)
shade(rep(mycisim[,2],each=2),rep(mycisim[,3],each=2),1:(2*(dim(myci)[1]))+0.5,1:(2*(dim(myci)[1]))+0.5,col="lightgrey",border=NA)
shade(rep(mycisim[,4],each=2),rep(mycisim[,5],each=2),1:(2*(dim(myci)[1]))+0.5,1:(2*(dim(myci)[1]))+0.5,col="lightgrey",border=NA)


for(ii in 1:dim(myci)[1]) {
  #axislab <- cbind(axislab,k)
  if(myci[ii,1]==1) {points(0,2*ii,col=color[ii],pch=8,cex=1)}
  if(mycisim[ii,1]==1){points(0,2*ii,col="lightgrey",pch=8,cex=1.5)}
  yval = rep(2*ii,length(myci[ii,])-1)
  xval = myci[ii,-1]
  if(length(xval)>2) {xval = c(xval[1:2],NA,xval[3:4]); yval = c(yval[1:2],NA,yval[3:4])}
  lines(xval,yval,col=color[ii],lty=linetype[ii],lwd=2)
  points(myest[ii],2*ii,col=color[ii],pch=20,cex=1)
}


axis(2,2*non0,names(non0),cex.axis=.6,las=T)


dev.off()




result.ci.cov0
