source("../Rscripts/extract.features.R")
source("../Rscripts/procrustes.R")
source("../Rscripts/transport.R")
source("../Rscripts/delay.reconstruction.R")

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


fly.ids <- c()
for(i in 1:length(flies)){
 source(flies[[i]])

  #extract features
  WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
    nRuns=5, std=0.1, lT=5, uT=-1)

  #extract segments from features
  Stmp <- extract.all.segments(WT, k=delay)
  fly.ids <- c(fly.ids, rep(i, length(Stmp)))
  Slist <- c(Slist, Stmp)

}
  
build.markov.transition <- function(clusters, time, nstates=max(clusters) ){

  M = matrix(0, nrow=nstates, ncol=nstates)
  for(i in 2:length(clusters)){
    if(time[i]-time[i-1] < 3){
      M[clusters[i-1], clusters[i] ]  = M[clusters[i-1], clusters[i] ] + 1
    }
  }

 for(i in 1:nrow(M)){
    M[i, ] = M[i, ] / sum(M[i, ])
  }
  M[is.na(M)]  = 0

#  M/length(clusters)
  M
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

  
nS = 6
nC = 6
sSeq = seq(-0.00001, 5.000001, length.out=6)
cSeq = seq(-pi, pi, length.out=6)


Mbefore = list()
Mduring = list() 
Mafter = list()

nCenters = (nS-1)*(nC-1)

for(i in 1:length(Slist)){
  Stmp  = Slist[[i]]
  cCuts = as.integer(cut(Stmp$curvatureMean, cSeq))
  sCuts = as.integer(cut(Stmp$lengths, sSeq))
  cluster = cCuts + ( (sCuts-1)*(nC-1) )

  index = Ilist[[i]]
  time = Stmp$C$timeToOdor
  o <- order(time)
  ind =  o[which(time < 0)] 
  Mbefore[[i]] = build.markov.transition(cluster[ind], time[ind], nCenters)
  ind = o[which(time >= 0 & time < 1000)]
  Mduring[[i]] = build.markov.transition(cluster[ind], time[ind], nCenters)
  ind = o[which( (time-min(time)) > 10800)] 
  Mafter[[i]] = build.markov.transition(cluster[ind], time[ind], nCenters)
}


meanBefore = list()
meanDuring = list()
meanAfter = list()
for(k in 1:length(fly.type)){
  ind = which(fly.ids == k)
  meanBefore[[k]] = matrix(0, nrow=nCenters, ncol=nCenters) 
  meanDuring[[k]] = matrix(0, nrow=nCenters, ncol=nCenters) 
  meanAfter[[k]] = matrix(0, nrow=nCenters, ncol=nCenters)
    
  for(i in ind){
    meanBefore[[k]] = meanBefore[[k]] + Mbefore[[i]]
    meanDuring[[k]] = meanDuring[[k]] + Mduring[[i]]
    meanAfter[[k]] = meanAfter[[k]] + Mafter[[i]]
  } 
  
  for(i in 1:nrow(meanBefore[[k]])){
    meanBefore[[k]][i, ] =  meanBefore[[k]][i, ] / sum( meanBefore[[k]][i, ] )
    meanDuring[[k]][i, ] =  meanDuring[[k]][i, ] / sum( meanDuring[[k]][i, ] )
    meanAfter[[k]][i, ] =  meanAfter[[k]][i, ] / sum( meanAfter[[k]][i, ] )
  }
  meanBefore[[k]][is.na(meanBefore[[k]])]  = 0
  meanDuring[[k]][is.na(meanDuring[[k]])]  = 0
  meanAfter[[k]][is.na(meanAfter[[k]])]  = 0

#meanBefore[[k]] = meanBefore[[k]]/length(ind)
# meanDuring[[k]] = meanDuring[[k]]/length(ind)
# meanAfter[[k]] = meanAfter[[k]]/length(ind)
}

layout( matrix(1:9, nrow=3))


library(RColorBrewer)
pal = brewer.pal(n=9, "Reds")
ramp=colorRamp(pal)
for(k in 1:length(fly.type)){
  image( z=t(meanBefore[[k]]), col=rgb(ramp( ((1:101)-1)/100)/255), asp=1 )
  image( z=t(meanDuring[[k]]), col=rgb(ramp( ((1:101)-1)/100)/255), asp=1 )
  image( z=t(meanAfter[[k]]), col=rgb(ramp( ((1:101)-1)/100)/255), asp=1 )
}


