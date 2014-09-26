source("../Rscripts/extract.features.R")
source("../Rscripts/procrustes.R")
source("../Rscripts/transport.R")
source("../Rscripts/delay.reconstruction.R")
source("../Rscripts/markov.R")

#to install packages uncomment:
#install.packages("shapes")
#install.packages("abind")
#install.packages("circular")

#load WT file names, adjust path in WT_ACV0-files.R tp point to the correct
#directoty. To convert .mat files (rim information files) to csv use
#rim-mat-2-csv.m in OctaveScripts (should be compatible with Matlab

flies = c("../Rscripts/WT_ACV0-files.R", "../Rscripts/Orco_ACV0-files.R", "../Rscripts/IR8a1_ACV0-files.R")
fly.type = c("WT", "Orco", "IR8a1")
delay=1


Slist <- list()

lT=4
fly.ids <- c()
for(i in 1:length(flies)){
 source(flies[[i]])

  #extract features
  WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
    nRuns=10, std=0.1, lT=lT, uT=-1)

  #extract segments from features
  Stmp <- extract.all.segments(WT, k=delay)
  fly.ids <- c(fly.ids, rep(i, length(Stmp)))
  Slist <- c(Slist, Stmp)

}

times <- c()
for(i in 1:length(Slist) ){
  times <- c(times, Slist[[i]]$C$timeToOdor)
}

maxIndex = 0
Ilist <-list()
for(k in 1:length(Slist) ){
  Ilist[[k]] = (maxIndex+1):(maxIndex+nrow(Slist[[i]]$C)) 
  maxIndex = maxIndex + nrow(Slist[[i]]$C) 

}

  
nS = 11
nC = 11
sSeq = seq(-0.00001, lT+ 0.000001, length.out=nS)
cSeq = seq(-pi, pi, length.out=nC)

curvature = rep( (cSeq[1:(nC-1)]+cSeq[2:nC])/2, nS-1)
speed = as.vector(t( matrix( (sSeq[1:(nS-1)]+sSeq[2:nS])/2,nrow=nC-1,
    ncol=nS-1)) )

Mbefore = list()
Mduring = list() 
Mafter = list()

nCenters = (nS-1)*(nC-1)

seq.length = 1000
all <- list()
index = 1
flytype = c()
flyid = c()
times = c();
for(i in 1:length(Slist)){
  Stmp  = Slist[[i]]
  cCuts = as.integer(cut(Stmp$curvatureMean, cSeq))
  sCuts = as.integer(cut(Stmp$lengths, sSeq))
  cluster = cCuts + ( (sCuts-1)*(nC-1) )

  intervals <- seq(1, length(cluster), by = seq.length)
  time = Stmp$C$timeToOdor
  o <- order(time)
  for(j in 2:length(intervals)){
    
    ind =  o[intervals[j-1]:intervals[j]] 
    times = c(times, mean(time[ind]))
    flyid = c(flyid, i)
    flytype = c(flytype, fly.ids[i])
    all[[index]] = build.markov.transition( cluster[ind], time[ind], nCenters )
    index = index +1
  }
}


X <- matrix(nrow = length(all), ncol=nCenters^2)
for(i in 1:length(all) ){
  X[i, ] = as.vector(all[[i]])
}


