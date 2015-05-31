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
source("../Rscripts/WT_ACV0-files.R")


#extract features
WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
  nRuns=10, std=0.1, lT=5, uT=-1)

#extract segments from features
Slist <- extract.all.segments(WT, k=100, offset=50)

flies = order(runif(30))[1:5]
Z <- list()
times <- list()
ids <- list()
odor <- list()
odor2 <- list()
index = 1
for( id in flies ){
  for(j in 1:nrow(Slist[[id]]$lengths) ){
    Z[[index]] = cbind( Slist[[id]]$curvatureMean[j, ], Slist[[id]]$lengths[j, ] )
    ids[[index]] = id
    times[[index]] = mean( Slist[[id]]$time[j, ] )
    odor[[index]] = Slist[[id]]$C$maxOdor[j]
    odor2[[index]] = odor[[index]]
    if(odor[[index]] > 1 & odor[[index]] < 5){
      odor[[index]] = 2
    }
    index = index + 1
  }
}

  
trp <- pairwise.transport(Xin = Z, eps=0.01, scale=-1, store.plan=F)




library(RColorBrewer)

a = c()
index = 0
k = length(flies)
g = unlist(ids)
np = c() 
for(id in flies ){
  n = sum(g==id)
  np = c(np, n)
  a = c(a, seq(index/k, (index+1)/k-0.001, length.out = n ) )
  index = index +1 
}

Ds = trp$D
cnp = c(0, cumsum(np))
for(i in 2:length(cnp)){
  xi = cnp[i-1]:cnp[i]
  for(j in 2:length(cnp)){
    yi = cnp[j-1]:cnp[j]
    Ds[xi, yi] = Ds[xi, yi] / max( Ds[xi, yi]) 
  }  
}

pal = brewer.pal("YlOrRd", n=9)
image(x=a, y=a, Ds, col=pal, bty="n", xlab="", ylab="", useRaster=F,  asp=1,
    xaxt="n", yaxt="n")

for(i in 2:length(odor)){
  if(odor[[i-1]] != odor[[i]]){
    x = (a[i] + a[i-1])/2
    lines(c(0, 1), c(x, x), col="#0000FF99")
    lines(c(x, x), c(0, 1), col="#0000FF99")
  }
}

for(i in 1:(k-1)){
  x = i/k - 0.001
  lines(c(0, 1), c(x, x), col="black", lwd=2)
  lines(c(x, x), c(0, 1), col="black", lwd=2)
}









cmds = cmdscale(trp$D, k=100, eig=T )

library(Rtsne)
tsne = Rtsne(X=cmds$points, pca=F, check_duplicates=F)


if(F){

pal = brewer.pal("Dark2", n=length(flies))
g = as.factor( unlist(ids) )

plot(cmds$points, col=pal[g], pch=19, bty="n", xlab="mds 1", ylab="mds 2", asp = 1) 
dev.copy2eps(file="../Rwasserstein/trp-seg-100-mds-12.eps")

plot(cmds$points[, c(2,3)], col=pal[g], pch=19, bty="n", xlab="mds 2", ylab="mds 3", asp = 1) 
dev.copy2eps(file="../Rwasserstein/trp-seg-100-mds-23.eps")

plot(cmds$points[, c(1,3)], col=pal[g], pch=19, bty="n", xlab="mds 1", ylab="mds 3", asp = 1) 
dev.copy2eps(file="../Rwasserstein/trp-seg-100-mds-13.eps")

plot(tsne$Y, col=pal[g], pch=19, bty="n", xlab="tsne 1", ylab="tsne  2", asp = 1) 
dev.copy2eps(file="../Rwasserstein/trp-seg-100-tsne.eps")
  
pal = c("darkblue", "darkred", "orange", "gold", "lightblue")
g = as.factor( unlist(odor2) )
  
plot(cmds$points, col=pal[g], pch=19, bty="n", xlab="mds 1", ylab="mds 2", asp = 1) 
dev.copy2eps(file="../Rwasserstein/trp-seg-100-all-odor-mds-12.eps")

plot(cmds$points[, c(2,3)], col=pal[g], pch=19, bty="n", xlab="mds 2", ylab="mds 3", asp = 1) 
dev.copy2eps(file="../Rwasserstein/trp-seg-100-odor-mds-23.eps")

plot(cmds$points[, c(1,3)], col=pal[g], pch=19, bty="n", xlab="mds 1", ylab="mds 3", asp = 1) 
dev.copy2eps(file="../Rwasserstein/trp-seg-100-odor-mds-13.eps")

plot(tsne$Y, col=pal[g], pch=19, bty="n", xlab="tsne 1", ylab="tsne  2", asp = 1) 
dev.copy2eps(file="../Rwasserstein/trp-seg-100-odor-tsne.eps")

}
