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

  
nS = 12
nC = 12
sSeq = seq(-0.00001, lT+ 0.000001, length.out=nS)
cSeq = seq(-pi, pi, length.out=nC)

curvature = rep( (cSeq[1:(nC-1)]+cSeq[2:nC])/2, nS-1)
speed = as.vector(t( matrix( (sSeq[1:(nS-1)]+sSeq[2:nS])/2,nrow=nC-1,
    ncol=nS-1)) )



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
  ind = o[which(time >= 0 & time < 2500)]
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

layout( t( matrix(1:3, ncol=3) ))


library(RColorBrewer)
pal = c("#FFFFFF", brewer.pal(n=9, "Reds")[4:9])
ramp=colorRamp(pal)
cols = rgb(ramp( ((1:200)-1)/199)/255)

palD = rev(brewer.pal(n=9, "RdBu"))[c(1,2,3,3,5,7,7,8,9)]
rampD=colorRamp(palD)
colsD = rgb(rampD( ((1:200)-1)/199)/255)

for(k in 1:length(fly.type)){
 if(F){
  off = 0.5/(nCenters+1)
  image( z=t(meanBefore[[k]]), col=cols, asp=1, zlim=c(0,1),
      bty="n", xaxt="n", yaxt="n")
 abline(v=seq(-off,1+off, length.out=(nS-1)*(nC-1)+1), col="lightgray", lty=3)
 abline(h=seq(-off,1+off, length.out=(nS-1)*(nC-1)+1), col="lightgray", lty=3)
 abline(v=seq(-off,1+off, length.out=nS), col="black", lty=1)
 abline(h=seq(-off,1+off, length.out=nS), col="black", lty=1)
 title(sprintf("Transitions %s Before Odor", fly.type[k] ))

  image( z=t(meanDuring[[k]]), col=cols, asp=1 , zlim=c(0,1), bty="n", xaxt="n",
      yaxt="n")
 abline(v=seq(-off,1+off, length.out=(nS-1)*(nC-1)+1), col="lightgray", lty=3)
 abline(h=seq(-off,1+off, length.out=(nS-1)*(nC-1)+1), col="lightgray", lty=3)
 abline(v=seq(-off,1+off, length.out=nS), col="black", lty=1)
 abline(h=seq(-off,1+off, length.out=nS), col="black", lty=1)
  
 title(sprintf("Transitions %s During Odor", fly.type[k] ))
 }
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
  

  sB = rep(0, length(vB))
  sD = rep(0, length(vB))
  sA = rep(0, length(vB))
  
  msB = rep(0, nS-1)
  msD = rep(0, nS-1)
  msA = rep(0, nS-1)
 
  for( i in 1:(nS-1) ){
    start = ((i-1)*(nC-1) ) +1;
    end = start+nC-2
    sB[start:end] = sum(vB[start:end])
    sD[start:end] = sum(vD[start:end])
    sA[start:end] = sum(vA[start:end])
    msB[i] = sum(vB[start:end])
    msD[i] = sum(vD[start:end])
    msA[i] = sum(vA[start:end])
  }

  mcB = rep(0, nC-1)
  mcD = rep(0, nC-1)
  mcA = rep(0, nC-1)
  for( i in 1:(nS-1) ){
    start = ((i-1)*(nC-1) ) +1;
    end = start+nC-2
    mcB = mcB + vB[start:end]
    mcD = mcD + vD[start:end]
    mcA = mcA + vA[start:end]
  }
  v = max(c(sB, sD, sA))

  cB = vB / sB * v
  cD = vD / sD * v
  cA = vA / sA * v

  v = v*2;

  if(F){
  plot(sB, ylim = c(0, v), pch=19, col="purple", xaxt="n", bty="n", ylab="P",
    cex=0.7)
  points(sD, ylim = c(0, v), pch=19, col="orange", cex=0.5)
  points(sA, ylim = c(0, v), pch=19, col="blue", cex=0.5)
  
  
  points(cB, ylim = c(0, v), pch=19, col="purple", cex=0.5)
  points(cA, ylim = c(0, v), pch=19, col="blue", cex=0.5)
  points(cD, ylim = c(0, v), pch=19, col="orange", cex=0.5)
  
  xSeq1 = seq(nC-0.5, (nS-1)*(nC-1), by=nC-1)
  
  abline( v=xSeq1, col="black" )
  text( x=xSeq1-0.3, y=0.9*v, labels=sprintf("speed %.2f",sSeq[2:(length(sSeq)-1)]), srt=90 )
  
  xSeq2 = seq( nC/2, (nS-1)*(nC-1), by=nC-1)
  abline( v=xSeq2, col="gray" )
  text( x=xSeq2-0.3, y=0.9*v, labels="curvature -", srt=90 )
  text( x=xSeq2+0.2, y=0.9*v, labels="curvature +", srt=90 )

  title(sprintf("Stationary distributions %s (orange= during, purple = before, blue=after)", fly.type[k] ))
  }

  if(F){
  v = max(c(msB, msD, msA))
  plot(msB, ylim = c(0, v), pch=19, col="purple", xaxt="n", bty="n", ylab="P",
    cex=0.7, type="l", lwd=2)
  lines(msD, ylim = c(0, v), pch=19, col="orange", cex=0.5, lwd=2)
  lines(msA, ylim = c(0, v), pch=19, col="blue", cex=0.5, lwd=2)
  
  v = max(c(mcB, mcD, mcA))
  plot(mcB, ylim = c(0, v), pch=19, col="purple", xaxt="n", bty="n", ylab="P",
    cex=0.7, type="l", lwd=2)
  lines(mcD, ylim = c(0, v), pch=19, col="orange", cex=0.5, lwd=2)
  lines(mcA, ylim = c(0, v), pch=19, col="blue", cex=0.5, lwd=2)
  }



  pB = generate.path( meanBefore[[k]], c(nS/2, nC/2), 5000, speed, curvature )
  pB = prcomp(pB)$x
  pD = generate.path( meanDuring[[k]], c(nS/2, nC/2), 5000, speed, curvature )
  pD = prcomp(pD)$x
  pA = generate.path( meanAfter[[k]], c(nS/2, nC/2), 5000, speed, curvature )
  pA = prcomp(pA)$x
  rx = range(c(pB[,1], pD[,1], pA[,1]))
  ry = range(c(pB[,2], pD[,2], pA[,2]))
  plot(pB, xlim=rx, ylim=ry, type="l", col="purple", lwd=2, asp=1)
  lines(pD, xlim=rx, ylim=ry, col="orange", lwd=2, asp=1)
  lines(pA, xlim=rx, ylim=ry, col="blue", lwd=2, asp=1)
}


