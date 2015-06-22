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
fly_type = "WT"
seg_length = 30
seg_offset = 15

dir_name = sprintf("../Rwasserstein/Individuals2-%s-seg-%d-off-%d", fly_type, seg_length, seg_offset )
dir.create( dir_name )

if(F){
#extract features
WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
  nRuns=10, std=0.1, lT=5, uT=-1)

#extract segments from features
Slist <- extract.all.segments(WT, k=seg_length, offset=seg_offset)
}

flies = c(4,6,15,17,24) 
for( id in flies ){
 
  Z <- list()
  times <- list()
  timeToOdor <- list()
  odor <- list()
  odor2 <- list()
  index = 1
  for(j in 1:nrow(Slist[[id]]$lengths) ){
    Z[[index]] = cbind( Slist[[id]]$curvatureMean[j, ], Slist[[id]]$lengths[j, ] )
    times[[index]] = mean( Slist[[id]]$time[j, ] )
    timeToOdor[[index]] = Slist[[id]]$C$timeToOdor[j]
    odor[[index]] = Slist[[id]]$C$maxOdor[j]
    odor2[[index]] = odor[[index]]
    if(odor[[index]] > 1 & odor[[index]] < 5){
      odor[[index]] = 2
    }
    index = index + 1
  }

  
  trp <- pairwise.transport(Xin = Z, eps=0.01, scale=-1, store.plan=F)


  library(RColorBrewer)



  pal = brewer.pal("YlOrRd", n=9)
  image(trp$D, col=pal, bty="n", xlab="", ylab="", useRaster=F,  asp=1,
     xaxt="n", yaxt="n")

  n = ncol(trp$D)
  for(i in 2:length(timeToOdor)){
    if(timeToOdor[[i-1]] < 0 & timeToOdor[[i]] > 0 ){
      x = (i - 0.5) / n
      lines(c(0, 1), c(x, x), col="Magenta")
      lines(c(x, x), c(0, 1), col="Magenta")
    }
  }

  dev.copy2eps( file = sprintf("%s/fly-%d-pairwise.eps", dir_name, id) ) 









  cmds = cmdscale(trp$D, k=100, eig=T )

  library(Rtsne)
  tsne = Rtsne(X=cmds$points, pca=F, check_duplicates=F)


  cols = rep( rgb( t(col2rgb("gold2")/255) ), nrow(cmds$points) )
  g = unlist(odor)
  cols[g>1] = rgb( t(col2rgb("orange")/255) )
  ramp = colorRamp(c("magenta", "darkblue"))
  g = unlist(timeToOdor)
  g = g/max(g);
  cols[g>=0] = rgb(ramp(g[g>=0])/255)
 
  plot(cmds$eig, col="black", bty="n", xlab="Dimension",
      ylab="Eigenvalue", type="l", lwd=4) 
  dev.copy2eps(file = sprintf("%s/fly-%d-mds-eigs.eps", dir_name, id) )
    
  plot(cmds$eig[1:10], col="black", bty="n", xlab="Dimension",
      ylab="Eigenvalue", type="l", lwd=4) 
  dev.copy2eps(file = sprintf("%s/fly-%d-mds-eigs-1-10.eps", dir_name, id) )
  
  plot(cmds$points, col=cols, bty="n", xlab="mds 1", ylab="mds 2", asp = 1) 
  dev.copy2eps(file = sprintf("%s/fly-%d-mds-12.eps", dir_name, id) )

  plot(cmds$points[, c(2,3)], col=cols, bty="n", xlab="mds 2", ylab="mds 3", asp = 1) 
  dev.copy2eps(file = sprintf("%s/fly-%d-mds-23.eps", dir_name, id) )

  plot(cmds$points[, c(1,3)], col=cols, bty="n", xlab="mds 1", ylab="mds 3", asp = 1) 
  dev.copy2eps(file = sprintf("%s/fly-%d-mds-13.eps", dir_name, id) )

  plot(tsne$Y, col=cols, bty="n", xlab="tsne 1", ylab="tsne  2", asp = 1) 
  dev.copy2eps(file = sprintf("%s/fly-%d-tsne.eps", dir_name, id) )

}
