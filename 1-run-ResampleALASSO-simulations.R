#Author: Jessica Minnier
#Main file to run simulations

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


source("functions.R")

## Sets up parameter values for simulation

##----- SET THESE PARAMETERS
paramnum = 1
outputwd = getwd()
nptb = 500 		#number of betastars
nrep = 1500/numsets 		#number of data sets in each file, 1500/75 = 20 sets each

##-----LOAD PARAMETER VALUES
load("paramlist.RData")
myparams = paramlist[[paramnum]] #paramnum from args in command line
mysig0 = myparams[[1]]
cor.X = myparams[[2]]
sampsize = tmpn = myparams[[3]]
beta.true = myparams[[4]]
pz = length(beta.true)
setnum = myparams[[5]]
numsets = myparams[[6]]

fn.out = paste0("alasso.samp",sampsize,"ptb",nptb,"pz",pz,"b",
                sprintf("%1d",beta.true[4]*10),"sig",mysig0,"corX",cor.X,".set",setnum)
fn.out = file.path(outputwd,fn.out)

print(myparams)

mybictype="min"



fn.file = try(read.table(fn.out))
if(class(fn.file)=="try-error") {fn.file=matrix(1,nrow=1,ncol=1)}
print(c("nrow",nrow(fn.file)))

if(nrow(fn.file)<nrep) {
  system(paste("rm",fn.out))
  print(fn.out)
  
  for(m in 1:nrep)
  {
    mydata = SIM.FUN(sampsize,rho=cor.X);          
    means = apply(mydata,2,mean)
    mydata.c = mydata - VTM(means,sampsize)
    tmpratio=1
    
    #obtain alasso estimate
    beta.hat = Est.ALASSO(mydata.c,bictype=mybictype)
    beta.ptb.all = NULL
    #perturbation
    for(i in 1:nptb)
    {
      tmpG = rexp(nrow(mydata))
      means2 = apply(mydata*tmpG,2,mean)
      mydata.c2 = (mydata - VTM(means2,sampsize))*sqrt(tmpG)
      beta.ptb = Est.ALASSO(mydata.c2,rtn="notall",
                            s0=beta.hat[[1]][1],bictype=mybictype) #don't return ptb'd Ahat
      beta.ptb.all = rbind(beta.ptb.all,beta.ptb)
    }
    print(m)
    write(round(c(beta.hat$beta.out,beta.hat$Ahat,beta.hat$sighat2,beta.ptb.all),10),
          ncol=100000,file=fn.out,append=T)
  }
  
  ########### go to compress file
}#if file is not full
