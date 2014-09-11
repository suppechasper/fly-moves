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
seg.len <- 20
flies <- c("../Rscripts/WT_ACV0-files.R", "../Rscripts/Orco_ACV0-files.R", "../Rscripts/IR8a1_ACV0-files.R")
fly.type <- c("WT", "Orco", "IR8a1")

Slist <- c()
fly.type.id <- c()
fly.ids <- c()
for(k in 1:length(flies)){

  source(flies[[k]])

  #extract features
  WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
    nRuns=10, lT=5, uT=0.1)

  #extract segments from features
  Stmp <- extract.all.segments(WT, k=seg.len)
 
  fly.ids <- c(fly.ids, rep(k, length(Stmp)))

  for(i in 1:length(Stmp)){
    fly.type.id <- c(fly.type.id, rep(k, nrow(Stmp[[i]]$C)))
  }
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



#do joint procrustes analysis
gpa.joint <- procrustes.analysis.joint(Slist)

#extract the pca transformed coordinates of the segments for each subject
Z <- procrustes.extract.rawscores.joint(gpa.joint, Slist, 1:10)



  
build.markov.transition <- function(clusters, nstates=max(clusters) ){

  M = matrix(0, nrow=nstates, ncol=nstates)
  for(i in 2:length(clusters)){
    M[clusters[i-1], clusters[i] ]  = M[clusters[i-1], clusters[i] ] + 1
  }

 for(i in 1:nrow(M)){
    M[i, ] = M[i, ] / sum(M[i, ])
  }
M[is.na(M)]  = 0

#  M/length(clusters)
M
}

nCenters = 25
#cluster = kmeans(gpa.joint$rawscores[,1], nCenters, nstart=10)

pc = 2
factor=3
s = seq(min(gpa.joint$rawscores[,pc])*(1.0001), max(gpa.joint$rawscores[,pc])*(1.0001), length.out =
      nCenters+1 )
cluster = as.integer( cut(gpa.joint$rawscores[,pc], breaks =s ))


Mbefore = list()
Mduring = list() 
Mafter = list()


for(i in 1:length(Slist)){
  Stmp  = Slist[[i]]
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


