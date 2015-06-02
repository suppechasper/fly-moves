source("../Rscripts/extract.features.R")
source("../Rscripts/procrustes.R")
source("../Rscripts/transport.R")
source("../Rscripts/delay.reconstruction.R")

#to install packages uncomment:
#install.packages("shapes")
#install.packages("abind")
#install.packages("circular")

#load WT file names, adjust path in WT_ACV0-files.R tp point to the correct
#directory. To convert .mat files (rim information files) to csv use
#rim-mat-2-csv.m in OctaveScripts (should be compatible with Matlab
source("../Rscripts/WT_ACV0-files.R")
fly_type = "WT"


#extract features
WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles, nRuns=10)

#extract segments from features
Slist <- extract.all.segments(WT, k=1)

Z = list()
for( i in 1:length(Slist) ){
  Z[[i]] =  scale( cbind(Slist[[i]]$curvatureMean/3.15, Slist[[i]]$lengths/5 ) ) 
}

D <- c()
Olist <- list()
paths <- list()

for(i in 1:length(Z)){
  #extract all sgements based on odor condition
  O <- extract.condition.odor(Z[[i]], Slist[[i]]$C) 

  #put all odors in a list
  Olist <- c(Olist, O)
}


trp <- pairwise.transport( Xin = Olist, stop=3, eps=0.4, scale=-1, store.plan=F,
    rFactor=1, p=1 )



dir_name = sprintf("../Rwasserstein/Odors-%s", fly_type)
dir.create( dir_name )

pal = brewer.pal("YlOrRd", n=9)

image(trp$D, col=pal, bty="n", xlab="", ylab="", useRaster=F,  asp=1,
    xaxt="n", yaxt="n")
a = 0.5/151
for(i in 1:29 ){
  x = i/30 *(1 + 1/151) - 0.5/151 

  lines(c(0-a, 1+a), c(x, x), col="black", lwd=1)
  lines(c(x, x), c(0-a, 1+a), col="black", lwd=1)
}
dev.copy2eps( file = sprintf("%s/pairwise.eps", dir_name) ) 


Ds = trp$D
for(i in seq(5, 150, by=5) ){
  xi = (i-5):i
  for(j in seq(5, 150, by=5) ){
    yi = (j-4):j
    Ds[xi, yi] = Ds[xi, yi] / max( Ds[xi, yi]) 
  }  
}

image(Ds, col=pal, bty="n", xlab="", ylab="", useRaster=F,  asp=1,
    xaxt="n", yaxt="n")
for(i in 1:29 ){
  x = i/30 *(1 + 1/151) - 0.5/151 
  lines(c(0-a, 1+a), c(x, x), col="black", lwd=1)
  lines(c(x, x), c(0-a, 1+a), col="black", lwd=1)
}
dev.copy2eps( file = sprintf("%s/pairwise-scaled.eps", dir_name) ) 





cmds = cmdscale(trp$D, k=100, eig=T )

library(Rtsne)
tsne = Rtsne(X=cmds$points, pca=F, check_duplicates=F)
 
#color by odor
ramp = colorRamp(c("magenta", "darkblue"))
colsAll = rep( c("gold2", rgb(ramp(seq(0,1, length.out=4) )/255)), 30)
cols = colsAll

addLines <- function(X){
  for(i in seq(2, 60, by=2) ){
    lines(X[c(i-1, i), ], col="#000000", lwd=1)
  }
}

plot(cmds$eig, col="black", bty="n", xlab="Dimension",
    ylab="Eigenvalue", type="l", lwd=4)
dev.copy2eps(file = sprintf("%s/mds-eigs.eps", dir_name) )
    
plot(cmds$eig[1:10], col="black", bty="n", xlab="Dimension",
      ylab="Eigenvalue", type="l", lwd=4) 
dev.copy2eps(file = sprintf("%s/mds-eigs-1-10.eps", dir_name) )


type = c("", "centered")
for(k in 1:length(type)){
  points= cmds$points
  pointsTsne = tsne$Y
  if( k > 1){
    for(i in seq(5, 150, by=5) ){
      points[(i-4):i, ] = points[(i-4):i, ] - 
          t( matrix(points[i-4, ], nrow=ncol(points), ncol=5 ))
      pointsTsne[(i-4):i, ] = pointsTsne[(i-4):i, ] - 
          t( matrix(pointsTsne[i-4, ], nrow=ncol(pointsTsne), ncol=5 ))
    }
  }

  ind = as.vector(rbind(1+5*(0:29), 2+5*(0:29)))

  points = points[ ind, ] 
  pointsTsne = pointsTsne[ ind, ]
  cols = colsAll[ind]

  plot(points, col=cols, bty="n", xlab="mds 1", ylab="mds 2", asp = 1) 
  addLines(points[, c(1,2)] )
  points(points, col=cols) 
  dev.copy2eps(file = sprintf("%s/mds-%s-12.eps", dir_name, type[k]) )

  plot(points[, c(2,3)], col=cols, bty="n", xlab="mds 2", ylab="mds 3", asp = 1) 
  addLines(points[, c(2,3)] )
  points(points[, c(2,3)], col=cols) 
  dev.copy2eps(file = sprintf("%s/mds-%s-23.eps", dir_name, type[k]) )

  plot(points[, c(1,3)], col=cols, bty="n", xlab="mds 1", ylab="mds 3", asp = 1) 
  addLines(points[, c(1,3)] )
  points(points[, c(1,3)], col=cols) 
  dev.copy2eps(file = sprintf("%s/mds-%s-13.eps", dir_name, type[k]) )


  plot(pointsTsne, col=cols, bty="n", xlab="tsne 1", ylab="tsne  2", asp = 1) 
  addLines(pointsTsne )
  points(pointsTsne, col=cols) 
  dev.copy2eps(file = sprintf("%s/tsne-%s.eps", dir_name, type[k]) )

}



