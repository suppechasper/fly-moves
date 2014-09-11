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

O.all <- list()
W.all <- list()
O.names <- c()

for(i in 1:length(flies)){
 source(flies[[i]])

  #extract features
  WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
    nRuns=10, std=0.1, lT=5, uT=-1)

  #extract segments from features
   Slist <- extract.all.segments(WT, k=delay)

  Z <- list()
  W <- list()
  for( j in 1:length(Slist) ){
    Z[[j]] = cbind(Slist[[j]]$curvatureMean, Slist[[j]]$lengths)
    W[[j]] <- rep(1/length(Z[[j]]), length(Z[[j]]))
  }

  #extracct all sgements based on odor condition
  O <- extract.condition.odor.all(Z, Slist) 
  names = names(O)

  O.all = c(O.all,  O)
  W.all = c(W.all,  W)
  O.names = c(O.names, sprintf("%s-%s", fly.type[i], names))


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





sSeq = seq(-0.00001, 5.000001, length.out=6)
cSeq = seq(-pi, pi, length.out=6)

cluster = as.integer( cut(gpa.joint$rawscores[,pc], breaks =s ))

Mbefore = list()
Mduring = list() 
Mafter = list()


for(i in 1:length(Slist)){
  Stmp  = Slist[[i]]
  cCuts = cut(Stmp$curvatureMean, cSeq)
  sCuts = cut(Stmp$lengths, sSeq)
  cluster = cCuts + ( (sCuts-1)*max(cCuts) + 1)
  nCenters = max(cluster)

  index = Ilist[[i]]
  time = Stmp$C$timeToOdor
  o <- order(time)
  ind = index[ o[which(time < 0)] ]
  Mbefore[[i]] = build.markov.transition(cluster[ind], nCenters)
  ind = index[ o[which(time >= 0 & time < 1000)] ]
  Mduring[[i]] = build.markov.transition(cluster[ind], nCenters)
  ind = index[ o[which( (time-min(time)) > 10800)] ]
  Mafter[[i]] = build.markov.transition(cluster[ind], nCenters)
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

layout( matrix(1:12, nrow=3))
  d = dim( gpa.joint$mshape )
  dir = factor* gpa.joint$pcasd[pc] *matrix(gpa.joint$pcar[,pc] , nrow=d[1])
  pc.segment.plot(gpa.joint$mshape, dir, factor, pc)


library(RColorBrewer)
pal = brewer.pal(n=9, "Reds")
ramp=colorRamp(pal)
x = (s[1:nCenters] + s[-1]) /2
for(k in 1:length(fly.type)){
  image( z=t(meanBefore[[k]]), x=x, y=x, col=rgb(ramp( ((1:101)-1)/100)/255), asp=1 )
  image( z=t(meanDuring[[k]]), x=x, y=x, col=rgb(ramp( ((1:101)-1)/100)/255), asp=1 )
  image( z=t(meanAfter[[k]]), x=x, y=x, col=rgb(ramp( ((1:101)-1)/100)/255), asp=1 )
}


