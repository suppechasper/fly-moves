##
#frame odor on:   5400
#frame odor off: 10801
#fly dimensions 25 x 12 pixels
#instantaneous velocity < 3 pixels /s , values abouve 4 p/s probably errorenous
#peak fluctation is 0.15 pixels, sd 0.09 pixels




conditional.procrustes <- function(S){
  library(shapes)

  C = as.data.frame( S$conditions )

  ind = which(C$odorOn == 1 & C$odorOff == 0 & C$innerRim==1 & C$middle==0 &
      C$outerRim ==0 & C$id == -1 & S$nSamples > 0)
  X <- extract.procrustes.segments.block(S$samples[ind])
  gpa <- procGPA(X, scale=F, reflect=T, proc.output=T, distances=F, pcaoutput=F)
  save(gpa, ind, file="segs-20-gpa-odor-inside.Rdata")

  ind = which( C$middle==1 & C$innerRim==0 & C$outerRim==0 & C$id==-1, S$nSamples > 0)
  X <- extract.procrustes.segments.block(S$samples[ind])
  gpa <- procGPA(X, scale=F, reflect=T, proc.output=T, distances=F, pcaoutput=F)
  save(gpa, ind, file="segs-20-gpa-middle.Rdata")

  ind = which( C$middle==0 & C$innerRim==0 & C$outerRim==1 & C$id==-1 & S$nSamples > 0)
  X <- extract.procrustes.segments.block(S$samples[ind])
  gpa <- procGPA(X, scale=F, reflect=T, proc.output=T, distances=F, pcaoutput=F)
  save(gpa, ind, file="segs-20-gpa-atOuterRim.Rdata")

}





extract.procrustes.segments.block <- function(data){
  library(shapes)
  library(abind)
  
  l = ncol(data[[1]])
  s = l/2
  segs <- c()
  for(seg in data){
     x = abind::abind(seg[,1:s], along=1)
     y = abind::abind(seg[,(s+1):l], along=1)
     z = aperm( abind::abind(x, y, along=3), c(2,3,1) )
     segs = abind::abind(segs, z, along=3)
  }
  segs
}





extract.procrustes.segments <- function(data){
  library(shapes)
  library(abind)
  k <- length( data[[1]]$segsX1 )
  segs <- c()
  for(seg in data){
     x = abind::abind(seg$segsX1, along=1)
     y = abind::abind(seg$segsX2, along=1)
     z = aperm( abind::abind(x, y, along=3), c(2,3,1) ) 
     segs = abind::abind(segs, z, along=3)
  }
  segs
}




