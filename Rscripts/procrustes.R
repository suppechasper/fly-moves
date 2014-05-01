##
#frame odor on:   5400
#frame odor off: 10801
#fly dimensions 25 x 12 pixels
#instantaneous velocity < 3 pixels /s , values abouve 4 p/s probably errorenous
#peak fluctation is 0.15 pixels, sd 0.09 pixels


procrustes.extract.rawscores.joint <-function(gpa, Slist){
  
  Z = list()
  index = 1
  for(i in 1:length(Slist)){
    n = nrow(Slist[[i]]$C)
    Z[[i]] = gpa$rawscores[index:(index+n-1), ]
    index = index + n
  }

  Z
} 

procrustes.extract.tangent.joint <-function(gpa, Slist){
  
  Z = list()
  index = 1
  for(i in 1:length(Slist)){
    n = nrow(Slist[[i]]$C)
    Z[[i]] = t(gpa$tan[,index:(index+n-1)])
    index = index + n
  }

  Z
} 

procrustes.extract.rotated.joint <-function(gpa, Slist){
  
  Z = list()
  index = 1
  for(i in 1:length(Slist)){
    n = nrow(Slist[[i]]$C)
    Z[[i]] = gpa$rotated[,,index:(index+n-1)]
    d = dim(Z[[i]])
    dim(Z[[i]]) = c(d[1]*d[2], d[3])
    Z[[i]] = t(Z[[i]])

    index = index + n
  }

  Z
}



extract.procrustes.block.from.segments <- function(Slist){
 
  Pall = c()
  for(S in Slist){
   P  = extract.procrustes.block.from.segment(S)
   Pall = abind::abind(Pall, P, along=3)
  }

  Pall
}



extract.procrustes.block.from.segment <- function(S){
  library(abind)
  X = abind::abind(S$x, S$y, along=3)
  aperm(X, c(2,3,1) ) 
}




procrustes.analysis <- function(S, scale=F){
  P = extract.procrustes.block.from.segment(S)
  procGPA(P, scale=scale, pcaoutput=T, proc.output=T, distances=F)
}




procrustes.analysis.joint <- function(Slist){
  P = extract.procrustes.block.from.segments(Slist)
  procGPA(P, scale=F, pcaoutput=T, proc.output=T)
}



procrustes.analysis.transport.plan <- function(trp){
  n = length(trp$cost)
  X = trp$from[[n]][trp$map[[n]][,1], ]  - trp$to[[n]][trp$map[[n]][,2], ]
#X = sweep(X, 1, trp$map[[n]][,3], "*")
  prcomp(X)  

   
}


#do procrustes per individual in list
procrustes.analysis.individual <- function(Slist){
  pga <- list()
  for(i in 1:length(Slist) ){
    pga[[i]] = procrustes.analysis(Slist[[i]])
  }
  pga
}




procrustes.align.individual <- function(pgas){

  P <- c()
  for( i in 1:length(pgas) ){
    P = abind::abind( P, pgas[[i]]$pcar, along=3) 

  }


  pga <- procGPA(P, scale=T, reflect=T, proc.output=T )
  
  mean <- matrix(0, nrow=nrow(pgas[[1]]$mshape), ncol=ncol(pgas[[1]]$mshape) )
  for( i in 1:length(pgas) ){
    mean = mean+pgas[[i]]$mshape
  }
  mean = mean/length( pgas )

  X = list()
  for( i in 1:length(pgas) ){
    X[[i]] = t(pga$rotated[,,i]) %*% pgas[[i]]$pcar %*% t(pgas[[i]]$rawscores) 
  }

  list(gpa=pga, mshape= mean, X=X )
  
}









procrustes.plot <- function(gpa, pc=1, factor=1){


  d = dim( gpa$mshape )
  dir = factor* gpa$pcasd[pc] *matrix(gpa$pcar[,pc] , nrow=d[1])
  
  pc.segment.plot(gpa$mshape, dir, factor, pc)

}

pc.segment.plot.vectorize <- function(pca, dim, factor, pc){
  mean = matrix(pca$center , nrow=dim)

  dir = factor* pca$sdev[pc] *matrix(pca$rotation[,pc] , nrow=dim)

  pc.segment.plot(mean, dir, factor, pc)
}

pc.segment.plot <- function(mean, dir, factor, pc){

  layout( matrix(1:3, nrow=1) )
  
  n = nrow(mean)
  xlim = range(mean[,1]) + range(dir[,1])
  ylim = range(mean[,2]) + range(dir[,2])
  plot(mean + dir, pch=19, xlim=xlim, ylim=ylim, asp=1, bty="n", xlab="",
      ylab="", cex.axis=2, cex=2, col="#00000050", type="l", lwd=5)
  points( (mean+dir)[2:n, ], pch=19, col="#00000090",cex=2)
  title(sprintf("Mean + %d sd * PC %d", factor, pc), cex.main=2)

  plot(mean, pch=19, xlim=xlim, ylim=ylim, asp=1, bty="n", xlab="",
      ylab="", cex.axis=2, cex=2, col="#00000050", type="l", lwd=5)
  points((mean)[2:n, ], pch=19, col="#00000090",cex=2)
  title("Mean", cex.main=2)
  
  plot(mean - dir, pch=19, xlim=xlim, ylim=ylim, asp=1, bty="n", xlab="",
      ylab="", cex.axis=2, cex=2, col="#00000050", type="l", lwd=5)
  points((mean - dir)[2:n, ], pch=19, col="#00000090",cex=2)
  title(sprintf("Mean - %d sd * PC %d", factor, pc), cex.main=2)
  
}





#does not make sense to do the analysis this way...
extract.invariant.segments <- function(Slist){
  Si <- c()
  for( i in 1:length(Slist) ){
    Si <- rbind(Si, cbind( t( apply( Slist[[i]]$rx, 1, cumsum)), 
                           t( apply( Slist[[i]]$ry, 1, cumsum))  ) )
  }

  Si
}


invariant.segments.pca.plot <- function(pca, dim=2, pc=1, factor=1){

  mean = matrix(pca$center, ncol=dim)
#mean = apply(mean, 2, cumsum)
  dir = factor*pca$sdev[pc] * matrix(pca$rotation[,pc], ncol=dim)
# dir = apply(dir, 2, cumsum)

  pc.segment.plot(mean, dir, sprintf("Mean + %d sd * PC %d", factor, pc))
}




