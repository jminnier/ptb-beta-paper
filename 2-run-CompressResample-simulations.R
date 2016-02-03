#Author: Jessica Minnier
#Compress all set files into one RData file
#Run after all 1-run-ResampleALASSo-simulations.R files have finished

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


##----- SET THESE PARAMETERS
paramnum = 1
outputwd = getwd()
fun.dir = outputwd
nptb = 500 		#number of betastars
yes.del = FALSE #delete RData file even if already created
setwd(outputwd)

##### LOAD PARAMETER VALUES
load("paramlist.RData")
load("paramlistsum.RData")
nptb = 500 									   	     #number of betastars
nrep = 1500/paramlist[[length(paramlist)]][[5]]	     #number of data sets in each file
numsets = paramlist[[length(paramlist)]][[5]]

myparams = paramlist.sum[[paramnum]] #paramnum from array number in lsf
mysig0 = myparams[[1]]
cor.X = myparams[[2]]
sampsize = tmpn = myparams[[3]]
beta.true = myparams[[4]]
pz = length(beta.true)

junk = NULL
tmp.fn = paste(outputwd,"alasso.samp",sampsize,"ptb",nptb,"pz",pz,"b",sprintf("%1d",beta.true[4]*10),"sig",mysig0,"corX",cor.X,sep="")
if(yes.del) {system(paste("rm ",tmp.fn,".RData",sep=""))}
for(setnum in 1:numsets) {
  
  tmp.fn.out = paste(tmp.fn,".set",setnum,sep="")
  if(!file.exists(paste(tmp.fn,".RData",sep=""))) {
    if(file.exists(tmp.fn.out)) {junk = rbind(junk,read.table(tmp.fn.out)  )}
    if( length(junk)>0) print(paste("setnum",setnum,"dim junk", dim(junk)[1],dim(junk)[2]))	
    if(setnum==numsets) {save(junk,file=paste(tmp.fn,".RData",sep=""))}
  }else{print(paste(tmp.fn,"exists."))}
}
