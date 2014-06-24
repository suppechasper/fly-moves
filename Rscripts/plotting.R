##
#frame odor on:   5400
#frame odor off: 10801
#fly dimensions 25 x 12 pixels
#instantaneous velocity < 3 pixels /s , values abouve 4 p/s probably errorenous
#peak fluctation is 0.15 pixels, sd 0.09 pixels

animate.cluster <- function(K, S, X, file="ac%.4d.png"){
  
  t = c()
  for(a in S$conds){
    t =c(t, mean(a[,3]))
  }


  for( i in 1:length(S$segs) ){
    png(sprintf(file, i), width=1500, height=1000)
   layout(mat=matrix(c(1,2,1,3,1,4), nrow=2, ncol=3))
#   plot(x=t, y=K, xlab="time", ylab="cluster id", pch=19, cex=0.5, col="#00000020")
#   abline(v=5400, col="gray")
#   abline(v=10800, col="gray")
    plot.cluster(K,S);

    abline(v=t[i], col="red")
    points(t[i], K[i], col="red", pch=19, cex=1.5)

    plot(X, t="l", asp=1, col="#00000020")
    lines( S$conds[[i]][, 4:5], col="red", lwd=3 )
    points( x=S$conds[[i]][1, 4], y =S$conds[[i]][1, 5], pch=19, cex=2)
    points( x=S$conds[[i]][nrow(S$conds[[i]]), 4], y
          =S$conds[[i]][nrow(S$conds[[i]]), 5], pch=6, cex=2)

    plot( S$conds[[i]][, 4:5], asp=1, type="l", col="red", lwd=3 )
    points( S$conds[[i]][, 4:5], asp=1, pch=19 )
    
    
    plot( S$segs[[i]][,1:2], asp=1, col="black", pch=19, xlim=c(0,6),
        ylim=c(-3.2,3.2) )
    points(x=mean(S$segs[[i]][,1]) , y = mean(S$segs[[i]][,2]), col="blue",
        pch=19, cex=2)
    lines(x=c(mean(S$segs[[i]][,1]),mean(S$segs[[i]][,1])) , y = c(
          mean(S$segs[[i]][,2]),sum(S$segs[[i]][,2])), col="red",lwd=2)
    abline(h=0)

    dev.off()

  }


} 



plot.cluster <- function(K, S){
  
  t = c()
  for(a in S$conds){
    t =c(t, mean(a[,3]))
  }

  cols = unlist( lapply(brewer.pal("Dark2", n=8), paste, "FF", sep="") )

  plot(x=t, y=K, xlab="time", ylab="cluster id", col=cols[K], pch=19, cex=0.5)
  abline(v=5400, col="gold")
  abline(v=10800, col="gold")
  t = c()
#for(i in 1:length(S$conds)){
#   a = S$conds[[i]]
#   lines(x =range(a[,3]), y=c(K[i], K[i]) )
#  }


  K.max = max(K)
  for(a in S$conds){
    y = rep(K.max/2, nrow(a))
    y[a[,1] >= -20 ] = 1.5
    y[a[,2] < 30 ] = K.max-0.5
    lines(a[,3], y, col="#00000020")
  }

  title("clusters")

} 


plot.cluster.density <- function(K, S, bw=200){
  library(RColorBrewer)

  t = c()
  for(a in S$conds){
    t =c(t, mean(a[,3]))
  }


  p = list()
  pMax = 0;
  for(i in unique(K)){
    tp = c();
    for(j in which(K==i) ){
      tp = c(tp, S$conds[[j]][,3])
    }
    p[[i]] = density(tp, n = 2^13, from=min(t), to=max(t), bw=bw,
        kernel="gaussian")

    pMax= max( max(p[[i]]$y), pMax)
  }

  plot(NULL, xlim=range(t), ylim=c(0, pMax), xlab="time", ylab="density", pch=19, cex=0.5)
  abline(v=5400, col="gold")
  abline(v=10800, col="gold")


  cols = unlist( lapply(brewer.pal("Dark2", n=8), paste, "90", sep="") )
#brewer.pal("Dark2", n=8)

  for(i in unique(K)){
    par(new=T)
    plot(p[[i]]$x, p[[i]]$y, col=cols[i], t="l", yaxt="n", xaxt="n", ylab="",
        xlab="", lwd=3);
  }
  
# par(new=T)
#  plot(NULL, xlim=range(t), ylim=c(0, 150), xaxt="n", yaxt="n", xlab="", ylab="", pch=19, cex=0.5)
#  for(a in S$conds){
#    lines(x=a[,3], y=a[,2], col="#00000020", t="l", yaxt="n", xaxt="n",
#        ylab="",xlab="")
#  }

  par(new=T)
   plot(NULL, xlim=range(t), ylim=c(0,1), yaxt="n", xaxt="n", ylab="",
        xlab="")
  for(a in S$conds){
    y = rep(0.5, nrow(a))
    y[a[,1] >= -10 ] = 0
    y[a[,2] < 20 ] = 1
    lines(a[,3], y, col="#00000020")
  }


  title("cluster density")
}





plot.clusters.distributions <- function(K, S){   

  library(MASS)

  n = ceiling(sqrt(max(K)))
  m =  ceiling(max(K)/n)
  layout( matrix(1:(m*n),m,n ))
  
  cols = unlist( lapply(brewer.pal("Dark2", n=8), paste, "30", sep="") )
  
  dK = vector("list", max(K))
  for( i in 1:length(S$segs) ){
    dK[[ K[i] ]]  = rbind(dK[[ K[i] ]], S$segs[[i]][,1:2])
  }

  for( i in 1:length(dK) ){
    plot( dK[[i]][,1:2], col=cols[i], pch=19, xlim=c(0,6),
        ylim=c(-3.2,3.2), xlab="", ylab="" )
      d1 = kde2d(x = dK[[i]][,1], y=dK[[i]][,2], n=100)
      contour(d1, levels=seq(0.01, max(d1$z), by=0.03), add=T, col="#00000050")
  }

}


plot.clusters.sample.segs <- function(K, S, nSegs=5){   

  library(MASS)

  n = ceiling(sqrt(max(K)))
  m =  ceiling(max(K)/n)
  layout( matrix(1:(m*n),m,n ))
  
  cols = unlist( lapply(brewer.pal("Dark2", n=8), paste, "90", sep="") )
  
  dK = vector("list", max(K))

  for( i in 1:length(S$segs) ){
    x = scale(S$conds[[i]][,4:5], center=T, scale=F)
    dK[[ K[i] ]]  = rbind( dK[[ K[i] ]], x )
  }

  sl = nrow(S$segs[[1]])
  for( i in 1:length(dK) ){
    l = nrow(dK[[i]]) / sl
    ind = floor(runif(nSegs)*l)
    plot( dK[[i]], col="#00000000", xlab="", ylab="", type="l", asp=1 )
    for(k in ind){
      q = k*sl+1
      lines( dK[[i]][q:(q+sl-1), ], col=cols[i] )
    }
  }

}



plot.curvature.length <- function(features){

  index1 = features$time < 5400
  index2 = features$time > 5400 & features$time < 10800
  index3 = features$time > 10800

   layout(mat=matrix(c(1,2,3), nrow=1, ncol=3))


plot(features$curvatureMean[index1],features$lengths[index1], pch=19,
    col="#00000010", xlab = "curvature", ylab = "speed")
d1 = kde2d(features$curvatureMean[index1],features$lengths[index1], n=100,
    h=0.6)
contour(d1, levels=seq(0.01, max(d1$z), by=0.03), add=T)
title("before")

plot(features$curvatureMean[index2],features$lengths[index2], pch=19,
    col="#00000010", xlab = "curvature", ylab = "speed")
d2 = kde2d(features$curvatureMean[index2],features$lengths[index2], n=100,
    h=0.6)
contour(d2, levels=seq(0.01, max(d2$z), by=0.03), add=T)
title("during")

plot(features$curvatureMean[index3],features$lengths[index3], pch=19,
    col="#00000010", xlab = "curvature", ylab = "speed")
d3 = kde2d(features$curvatureMean[index3],features$lengths[index3], n=100,
    h=0.6)
contour(d3, levels=seq(0.01, max(d3$z), by=0.03) , add=T)
title("after")

}


plot.radial.density <- function(features){
    index1 = features$time < 5400
  index2 = features$time > 5400 & features$time < 10800
  index3 = features$time > 10800

  d1 = density(x=features$dirim[index1], cut=0)
  d2 = density(x=features$dirim[index2], cut=0)
  d3 = density(x=features$dirim[index3], cut=0)
  d1$y = d1$y/(100-d1$x) *sum(1/(100-d1$x) )
  d2$y = d2$y/(100-d2$x) *sum(1/(100-d2$x) )
  d3$y = d3$y/(100-d3$x) *sum(1/(100-d3$x) )
1

  plot( d1$x, d1$y, type="l", col="black", ylim=c(0, max(c(d1$y, d2$y, d3$y))),
      xlab = "distance from inner rim", ylab="density" )
  lines(d2$x, d2$y, col="red")
  lines(d3$x, d3$y, col="blue")

  title("Normalized radial density")
  

}



plot.odor <- function(X, C, type="p"){
 cols = rep("#000000", nrow(C))
 ramp = colorRamp(c("red", "blue"))
 t = C$timeToOdor/max(C$timeToOdor)
 cols[t >= 0] =rgb( ramp(t[t>=0]^(1/2) )/255 ) 
 
 ramp2 = colorRamp(c("red", "black"))
 ind = C$timeToOdor > -60 & C$timeToOdor <0
 t = - C$timeToOdor[ind]
 t= t/max(t);

 
 plot(X, col=cols, pch=19, type=type)
}



plot3d.odor <- function(X, C, type="s", alpha=0.3){
 library(rgl)
 cols = rep("#000000", nrow(C))
 ramp = colorRamp(c("red", "darkblue"))
 t = C$timeToOdor/max(C$timeToOdor)
 cols[t >= 0] =rgb( ramp(t[t>=0]^(1/2) )/255 ) 
 
 ramp2 = colorRamp(c("red", "black"))
 ind = C$timeToOdor > -60 & C$timeToOdor <0
 t = - C$timeToOdor[ind]
 t= t/max(t);
 
 cols[ind] = rgb(ramp2(t)/255)


 plot3d(X, col=cols, type=type, r=1, asp="iso", alpha=alpha, lwd=2)

}



plot.location <- function(X, C){
 cols = c("black", "lightblue", "orange")
 plot(X, col=cols[C$maxPos], pch=19)
}




plot3d.location <- function(X, C){
 library(rgl)
 cols = c("black", "lightblue", "orange")
 plot3d(X, col=cols[C$maxPos], pch=19, type="s", r=0.025, asp="iso")





