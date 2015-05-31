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



pc.segment.plot <- function(mean, dir, factor, pc, cex.points=1){

  
  n = nrow(mean)
  xlim = range(mean[,1]) + range(dir[,1])
  ylim = range(mean[,2]) + range(dir[,2])

  plot(mean - dir, pch=19, xlim=xlim, ylim=ylim, asp=1, bty="n", xlab="",
      ylab="", cex.axis=2, cex=2, col="#00000050", type="l", lwd=5)
  points((mean - dir)[2:n, ], pch=19, col="#00000090",cex=cex.points)
  title(sprintf("Mean - %d sd * PC %d", factor, pc), cex.main=2)
  
  plot(mean, pch=19, xlim=xlim, ylim=ylim, asp=1, bty="n", xlab="",
      ylab="", cex.axis=2, cex=2, col="#00000050", type="l", lwd=5)
  points((mean)[2:n, ], pch=19, col="#00000090",cex=cex.points)
  title("Mean", cex.main=2)
  
  plot(mean + dir, pch=19, xlim=xlim, ylim=ylim, asp=1, bty="n", xlab="",
      ylab="", cex.axis=2, cex=2, col="#00000050", type="l", lwd=5)
  points( (mean+dir)[2:n, ], pch=19, col="#00000090",cex=cex.points)
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




#plot principal segments from procrustes analysis
procrustes.time.plot <- function(gpa, pc=1, factor=1, Z, times, ncuts=501,
    col="#000000", bw=100, add=F){


  library(KernSmooth)


  ctimes = cut(times, ncuts)
  levs = levels(ctimes)
  me = rep(0, length(levs))
  va = rep(0, length(levs))
  for(i in 1:length(levs) ){
    ind = which(ctimes == levs[i])
    me[i] = mean(Z[ind,pc])
    va[i] = var(Z[ind,pc])
  }


  x = seq(min(times), max(times), length.out=ncuts+1)
  x = (x[1:ncuts]+x[2:(ncuts+1)])/2
 
  rx = range(times)
  pme = locpoly(x=x, y=me, bandwidth=bw) 
  pva = locpoly(x=x, y=va, bandwidth=bw) 
  vam = 2*sqrt(max(pva$y))


  if(add){
   lines( pme$x, sqrt(pva$y), col=col, lwd=4 ) 
  }
  else{
    plot( pme$x, sqrt(pva$y), col=col, lwd=4, bty="n", type="l",
      ylim=c(0.5, 2)*range(sqrt(pva$y)) ) 
  }

#if(add){
#   lines( pme$x, pme$y+sqrt(pva$y), col=col, lwd=2 ) 
#  }
#  else{
#    plot( pme$x, pme$y+sqrt(pva$y), col=col, lwd=2, bty="n", type="l",
#      ylim=range(pme$y)+c(-vam, +vam) ) 
#  }
# lines(pme$x, pme$y         , col=col, lwd=5 ) 
#lines(pme$x, pme$y-sqrt(pva$y), col=col, lwd=2 ) 
#lines(ts, lp$fit + lp$fit.se, col=col, lwd=2) 
# lines(ts, lp$fit - lp$fit.se, col=col, lwd=2) 

}



#plot principal segments from procrustes analysis
procrustes.density.plot <- function(gpa, pc=1, factor=1, O,
    col=rep("#000000", length(O) ), n=1000){

  library(KernSmooth)

  par(mar=c(5,5,5,5))
  d = dim( gpa$mshape )
  dir = factor* gpa$pcasd[pc] *matrix(gpa$pcar[,pc] , nrow=d[1])
  
  layout( t(matrix(1:9, nrow=3)) )
  pc.segment.plot(gpa$mshape, dir, factor, pc)

  xlim <- c(0, 0)
  for(k in 1:length(O)){
    for(i in 1:length(O[[k]]) ){
       r <- range(O[[k]][[i]][,pc])
       xlim[1] <-  min(xlim[1], r[1])
       xlim[2] <-  max(xlim[2], r[2])
    }
  }


  ylim <- c(0,0)
  d <- list()
  means <- list()

  for(k in 1:length(O)){
    d[[k]] <- list()
    means[[k]] <- list()
    for(i in 1:length(O[[k]]) ){
     d[[k]][[i]] <- bkde(x=O[[k]][[i]][, pc], range.x=xlim)
      means[[k]][[i]] <- mean(x=O[[k]][[i]][,pc])
     ylim[2] = max(ylim[2], d[[k]][[i]]$y ) 
    }
  }


  for(i in 1:length(d[[1]]) ){
    yaxt = "n"
    ylab = ""
    if(i%%3==1){
      yaxt="s"
       ylab ="p"
    }
    plot( d[[1]][[i]], type="l", bty="n", lwd=2, col= col[1], xlim = xlim,
        ylim=ylim, cex.lab=1.5, xlab=names(O[[1]])[i], ylab=ylab, cex.axis=1.5,
        yaxt=yaxt)
    for(k in 2:length(O)){
      lines( d[[k]][[i]], lwd=2, col=col[k])
    }
    for(k in 1:length(O)){
      abline( v=means[[k]][[i]], lwd=1, col=col[k])
    }
    grid(col="#00000020", lty=1)
    abline(v=factor* gpa$pcasd[pc], col="#000000") 
    abline(v=-factor* gpa$pcasd[pc], col="#000000")

    if(i==1){
      legend(x="topright", col=cols, bty="n", legend=fly.type, cex=1.5, lwd=2,
        box.col="gray", bg="white")
    }
  }


}





procrustes.plot.shape <- function(gpa, pcs){
  shape = gpa$pcar[,1:length(pcs)] %*% pcs

  n = length(shape)
  plot(x=shape[1:(n/2)]+gpa$mshape[,1], y=shape[(n/2+1):n]+gpa$mshape[,2], type="l")
}







procrustes.plot.pcs <- function( gpa, index, O){

  layout(t(matrix(1:(length(index)*4), nrow=4)) )

  

  library(RColorBrewer)
  pal = brewer.pal(name="Dark2", n=length(O))

  for(i in 1:length(index) ){  
    pc = index[i]
    f=3
    if(pc==1){
      f=1;
    }
    d = dim( gpa$mshape )
    dir = f * gpa$pcasd[pc] *matrix(gpa$pcar[,pc] , nrow=d[1])
    pc.segment.plot(gpa$mshape, dir, f, pc)

    #plot variations
    se = c()
    counts = c()
    nperiods = 0
    for(i in 1:length(O)){
      nperiods = length(O[[i]])
        for(k in 1:length(O[[i]]) ){
          se = c(se, sum(O[[i]][[k]][,pc]^2))
            counts = c(counts, nrow(O[[i]][[k]]) )
        }
    } 
    Vtotal = sum(se)/(sum(counts)-1)
    
    Vperfly = c()
    for(i in seq(1, length(se), by=nperiods)){
      Vperfly = c(Vperfly,
      sum(se[i:(i+nperiods-1)])/(sum(counts[i:(i+nperiods-1)])-1) )
      i
    }
    Vperflyperodor = se/(counts-1)

    cols = rgb(0,0,0)
    cols = c(cols, pal)
    for(i in 1:length(O)){
      cols = c(cols, rep(pal[i], nperiods))
    }
    x = c(Vtotal, Vperfly, Vperflyperodor)
    names(x) = NULL
    c("total", "WT", "Orco", "IR8a", names(O[[1]]), names(O[[2]]),
          names(O[[3]]) ) 

    bp = barplot( x, col=cols)
    text(mp, par("usr")[3], labels = names(x), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
  }

}






procrustes.density.plot2 <- function(O, pc, factor, xlab=sprintf("PC %d", pc), legend=rep("", length(O)),
    col=rep("#000000", length(O) ), lty = rep(1, length(O)), n=1000, lwd=2){

  library(KernSmooth)
 
  xlim <- c(0, 0)
  for(k in 1:length(O)){
     r <- range(O[[k]][,pc])
     xlim[1] <-  min(xlim[1], r[1])
     xlim[2] <-  max(xlim[2], r[2])
  }


  ylim <- c(0,0)
  d <- list()
  means <- list()

  for(k in 1:length(O)){
    d[[k]] <- bkde(x=O[[k]][, pc], range.x=xlim)
    means[[k]] <- mean(x=O[[k]][,pc])
    ylim[2] = max(ylim[2], d[[k]]$y ) 
  }

  plot( NULL,  bty="n", xlim = xlim, ylim=ylim, xlab=xlab, ylab="density",
      cex.axis=1.5, cex.lab=1.5, yaxt="s")
  
  grid(col="#00000020", lty=1)
  for(i in 1:length(d) ){
    lines( d[[i]], type="l", lwd=lwd, col= col[i], lty=lty[i] )
    
    abline( v=means[[i]], lwd=1, col=col[i], lty=lty[i])
  }
  
  legend(x="topright", col=col, bty="n", legend=legend, cex=1.5, lwd=2,
    box.col="gray", bg="white", lty=lty)


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



