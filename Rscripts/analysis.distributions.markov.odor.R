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

lT=4
fly.ids <- c()
for(i in 1:length(flies)){
 source(flies[[i]])

  #extract features
  WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
    nRuns=10, std=0.1, lT=lt, uT=-1)

  #extract segments from features
  Stmp <- extract.all.segments(WT, k=delay)
  fly.ids <- c(fly.ids, rep(i, length(Stmp)))
  Slist <- c(Slist, Stmp)

}
  
build.markov.transition <- function(clusters, time, nstates=max(clusters) ){

  M = matrix(0, nrow=nstates, ncol=nstates)
  for(i in 2:length(clusters)){
    if(time[i]-time[i-1] < 2){
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
sSeq = seq(-0.00001, lT+ 0.000001, length.out=6)
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
pal = c("#FFFFFF", brewer.pal(n=9, "Reds")[2:9])
ramp=colorRamp(pal)
cols = rgb(ramp( ((1:200)-1)/199)/255)

palD = rev(brewer.pal(n=9, "RdBu"))[c(1,2,3,3,5,7,7,8,9)]
rampD=colorRamp(palD)
colsD = rgb(rampD( ((1:200)-1)/199)/255)

for(k in 1:length(fly.type)){
  off = 0.5/(nCenters+1)
  image( z=t(meanBefore[[k]]), col=cols, asp=1, zlim=c(0,1),
      bty="n", xaxt="n", yaxt="n")
 abline(v=seq(-off,1+off, length.out=(nS-1)*(nC-1)+1), col="gray", lty=3)
 abline(h=seq(-off,1+off, length.out=(nS-1)*(nC-1)+1), col="gray", lty=3)
 abline(v=seq(-off,1+off, length.out=nS), col="black", lty=1)
 abline(h=seq(-off,1+off, length.out=nS), col="black", lty=1)
 title(sprintf("Transitions %s Before Odor", fly.type[k] ))

  image( z=t(meanDuring[[k]]), col=cols, asp=1 , zlim=c(0,1), bty="n", xaxt="n",
      yaxt="n")
 abline(v=seq(-off,1+off, length.out=(nS-1)*(nC-1)+1), col="gray", lty=3)
 abline(h=seq(-off,1+off, length.out=(nS-1)*(nC-1)+1), col="gray", lty=3)
 abline(v=seq(-off,1+off, length.out=nS), col="black", lty=1)
 abline(h=seq(-off,1+off, length.out=nS), col="black", lty=1)
  
 title(sprintf("Transitions %s During Odor", fly.type[k] ))
  
#      yaxt="n")
#image( z=t(meanBefore[[k]] - meanDuring[[k]]), col=colsD, asp=1 ,
#zlim=c(-1,1), bty="n", xaxt="n", yaxt="n")
# abline(v=seq(-off,1+off, length.out=(nS-1)*(nC-1)+1), col="gray", lty=3)
# abline(h=seq(-off,1+off, length.out=(nS-1)*(nC-1)+1), col="gray", lty=3)
# abline(v=seq(-off,1+off, length.out=nS), col="black", lty=1)
# abline(h=seq(-off,1+off, length.out=nS), col="black", lty=1)
#  image( z=t(meanAfter[[k]]), col=cols, asp=1 , zlim=c(0,1), bty="n", xaxt="n",
  
  eB = eigen(t(meanBefore[[k]]))
  eD = eigen(t(meanDuring[[k]]))
  eA = eigen(t(meanAfter[[k]]))

  vB = abs( Re(eB$vector[, 1]) )
  vB = vB / sum(vB)
  vD = abs( Re(eD$vector[, 1]) )
  vD = vD / sum(vD)
  vA = abs( Re(eA$vector[, 1]) )
  vA = vA / sum(vA)
  v = max(c(vB, vD))*2

  plot(vB, ylim = c(0, v), pch=19, col="orange", xaxt="n", bty="n", ylab="P")
  points(vD, ylim = c(0, v), pch=19, col="purple")
  points(vA, ylim = c(0, v), pch=19, col="blue")
  xSeq1 = seq(nC-0.5, (nS-1)*(nC-1), by=nC-1)
  
  abline( v=xSeq1, col="black" )
  text( x=xSeq1-0.3, y=0.9*v, labels=sprintf("speed %.2f",sSeq[2:(length(sSeq)-1)]), srt=90 )
  
  xSeq2 = seq( nC/2, (nS-1)*(nC-1), by=nC-1)
  abline( v=xSeq2, col="gray" )
  text( x=xSeq2-0.3, y=0.9*v, labels="curvature -", srt=90 )
  text( x=xSeq2+0.2, y=0.9*v, labels="curvature +", srt=90 )

  title(sprintf("Stationary distributions %s (orange= before, purple = during, blue=after)", fly.type[k] ))
   
}


