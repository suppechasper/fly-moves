#---- Procrustes analysis of segments


#do procrustes analysis on the segments in S
procrustes.analysis <- function(S, scale=FALSE){
  library(shapes)
  P = extract.procrustes.block.from.segment(S)
  procGPA(P, scale=scale, pcaoutput=TRUE, proc.output=TRUE, distances=FALSE)
}



#do joint procrustes analysis on the list of segments in Slist
procrustes.analysis.joint <- function(Slist, scale=FALSE){
  library(shapes)
  P = extract.procrustes.block.from.segments(Slist)
  procGPA(P, scale=scale, pcaoutput=TRUE, proc.output=TRUE, distances=FALSE)
}


#do procrustes per individual in list
procrustes.analysis.individual <- function(Slist){
  library(shapes)
  pga <- list()
  for(i in 1:length(Slist) ){
    pga[[i]] = procrustes.analysis(Slist[[i]])
  }
  pga
}


#given a list fo segments build a the data structure for jointly doing a
#a procrustes analysis on all segments
extract.procrustes.block.from.segments <- function(Slist){
  library(abind)
  Pall = c()
  for(S in Slist){
   P  = extract.procrustes.block.from.segment(S)
   Pall = abind::abind(Pall, P, along=3)
  }

  Pall
}


#extract the data structure for procrustes analysis from the sgemenets in S
extract.procrustes.block.from.segment <- function(S){
  library(abind)
  X = abind::abind(S$x, S$y, along=3)
  aperm(X, c(2,3,1) ) 
}



#Extract per subject gpa results
procrustes.extract.rawscores.joint <-function(gpa, Slist, coords =
    1:nol(gpa$rawscores)){
  
  Z = list()
  index = 1
  for(i in 1:length(Slist)){
    n = nrow(Slist[[i]]$C)
    Z[[i]] = gpa$rawscores[index:(index+n-1), coords]
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

procrustes.extract.stdscores.joint <-function(gpa, Slist){  
  
  Z = list()
  index = 1
  for(i in 1:length(Slist)){
    n = nrow(Slist[[i]]$C)
    Z[[i]] = gpa$stdscores[index:(index+n-1), ]
    index = index + n
  }

  Z

}






#plot principal segments from procrustes analysis
procrustes.plot <- function(gpa, pc=1, factor=1){


  d = dim( gpa$mshape )
  dir = factor* gpa$pcasd[pc] *matrix(gpa$pcar[,pc] , nrow=d[1])
  
  layout( matrix(1:3, nrow=1) )
  pc.segment.plot(gpa$mshape, dir, factor, pc)

}

pc.segment.plot.vectorize <- function(pca, dim, factor, pc){
  mean = matrix(pca$center , nrow=dim)

  dir = factor* pca$sdev[pc] *matrix(pca$rotation[,pc] , nrow=dim)

  layout( matrix(1:3, nrow=1) )
  pc.segment.plot(mean, dir, factor, pc)
}



pc.segment.plot <- function(mean, dir, factor, pc){

  
  n = nrow(mean)
  xlim = range(mean[,1]) + range(dir[,1])
  ylim = range(mean[,2]) + range(dir[,2])

  plot(mean - dir, pch=19, xlim=xlim, ylim=ylim, asp=1, bty="n", xlab="",
      ylab="", cex.axis=2, cex=2, col="#00000050", type="l", lwd=5)
  points((mean - dir)[2:n, ], pch=19, col="#00000090",cex=2)
  title(sprintf("Mean - %d sd * PC %d", factor, pc), cex.main=2)
  
  plot(mean, pch=19, xlim=xlim, ylim=ylim, asp=1, bty="n", xlab="",
      ylab="", cex.axis=2, cex=2, col="#00000050", type="l", lwd=5)
  points((mean)[2:n, ], pch=19, col="#00000090",cex=2)
  title("Mean", cex.main=2)
  
  plot(mean + dir, pch=19, xlim=xlim, ylim=ylim, asp=1, bty="n", xlab="",
      ylab="", cex.axis=2, cex=2, col="#00000050", type="l", lwd=5)
  points( (mean+dir)[2:n, ], pch=19, col="#00000090",cex=2)
  title(sprintf("Mean + %d sd * PC %d", factor, pc), cex.main=2)
}


#plot principal segments from procrustes analysis
procrustes.qq.plot <- function(gpa, pc=1, factor=1, O, n=1000){


  par(mar=c(5,5,5,5))
  d = dim( gpa$mshape )
  dir = factor* gpa$pcasd[pc] *matrix(gpa$pcar[,pc] , nrow=d[1])
  
  layout( t(matrix(1:9, nrow=3)) )
  pc.segment.plot(gpa$mshape, dir, factor, pc)


  b1 = sort(O$Xbefore1[, pc])
  b1 = b1[seq(1, length(b1), length.out=n)]
  b2 = sort(O$Xbefore2[, pc])
  b2 = b2[seq(1, length(b2), length.out=n)]
  
  d1 = sort(O$Xduring1[, pc])
  d1 = d1[seq(1, length(d1), length.out=n)]
  d2 = sort(O$Xduring2[, pc])
  d2 = d2[seq(1, length(d2), length.out=n)]
  
  a1 = sort(O$Xafter1[, pc])
  a1 = a1[seq(1, length(a1), length.out=n)]
  a2 = sort(O$Xafter2[, pc])
  a2 = a2[seq(1, length(a2), length.out=n)]

  plot(b1, b2, xlab="Before 1", ylab="Before 2",
      pch=19, col="#00000030", cex.axis=2, cex.lab=2)
  abline(0,1, col="red")
  
  qqplot(b1, d1, xlab="Before 1", ylab="During 1",
      pch=19, col="#00000030", cex.axis=2, cex.lab=2)
  abline(0,1, col="red")
  
  qqplot(b1, d2, xlab="Before 1", ylab="During 2",
      pch=19, col="#00000030", cex.axis=2, cex.lab=2)
  abline(0,1, col="red")
  
  qqplot(b1, a1, xlab="Before 1", ylab="After 1",
      pch=19, col="#00000030", cex.axis=2, cex.lab=2)
  abline(0,1, col="red")
  
  qqplot(b1, a2, xlab="Before 1", ylab="After 2",
      pch=19, col="#00000030", cex.axis=2, cex.lab=2)
  abline(0,1, col="red")

}














#--- experimental ---


procrustes.difference.transport.plan <- function(trp){
  n = length(trp$cost)
  X = trp$from[[n]][trp$map[[n]][,1], ]  - trp$to[[n]][trp$map[[n]][,2], ]
#Xw = sweep(X, 1, trp$map[[n]][,3], "*")

  X
   
}


procrustes.analysis.transport.plan <- function(trp){
  n = length(trp$cost)
  X = trp$from[[n]][trp$map[[n]][,1], ]  - trp$to[[n]][trp$map[[n]][,2], ]
#X = sweep(X, 1, trp$map[[n]][,3], "*")
  prcomp(X)  

   
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




