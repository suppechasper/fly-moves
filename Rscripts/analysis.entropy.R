source("../Rscripts/extract.features.R")
source("../Rscripts/procrustes.R")
source("../Rscripts/transport.R")
source("../Rscripts/delay.reconstruction.R")


library(mmsentropy)
library(RColorBrewer)

#load WT file names, adjust path in WT_ACV0-files.R tp point to the correct
#directory. To convert .mat files (rim information files) to csv use
#rim-mat-2-csv.m in OctaveScripts (should be compatible with Matlab
flies = c( "../Rscripts/WT_ACV0-files.R", "../Rscripts/Orco_ACV0-files.R",  "../Rscripts/IR8a1_ACV0-files.R")
fly.type = c("WT", "Orco", "IR8a1")

scales = 20

#extract features

Elist <- list()
for(k in 1:length(flies)){
 source(flies[[k]])

  Flist <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles, nRuns=10)

  Slist <- extract.all.segments(Flist, k=1)

  Z <- list()
  for( i in 1:length(Slist) ){
    Z[[i]] = matrix( Slist[[i]]$lengths, ncol=1) 
#Z[[i]] = matrix(abs(Slist[[i]]$curvatureMean), ncol=1) 
  }

  O <- list()
  for(i in 1:length(Z) ){
    O[[i]] <- extract.condition.odor2( Z[[i]], Slist[[i]]$C ) 
  }


  E <- list()
  for(i in 1:length(O[[1]]) ){
    E[[i]] <- matrix(0, ncol=scales, nrow=length(O) )
    for(j in 1:length(O) ){
      if(length(O[[j]][[i]]) < 1000){
        E[[i]][j, ] = NA
      }
      else{
        E[[i]][j, ] <- multiscale.entropy(matrix(O[[j]][[i]], ncol=1), 2, scales,
           0.2)[1, ]
      } 
    }
  }

  Elist[[k]] = E
}


plot(NA, xlim=c(1, scales), ylim= c(1.2, 2.1), bty="n", xlab="scale", ylab="sample entropy" )  
cols = brewer.pal(n=3, "Dark2")
legend=c()
lty=c()
for(k in 1:length(Elist) ){
  E = Elist[[k]]

  means <- list()
  sdevs <- list()
  for(i in 1:length(E) ){
    means[[i]] = colMeans(E[[i]], na.rm=TRUE)
    sdevs[[i]] = apply(E[[i]],2, sd, na.rm=TRUE)
  }



  for(i in 1:length(means) ){
    lines(means[[i]], lwd=4, col=cols[i], lty=k)
    lines(means[[i]]+sdevs[[i]], lwd=1, col=cols[i], lty=k)
    lines(means[[i]]-sdevs[[i]], lwd=2, col=cols[i], lty=k)
  }

  #plot difference
  if(FALSE){
  means <- list()
  sdevs <- list()
  for(i in 1:length(E) ){
    means[[i]] = colMeans(E[[i]][, 2:scales] - E[[i]][, 1:(scales-1)])
    sdevs[[i]] = apply(E[[i]][, 2:scales] - E[[i]][, 1:(scales-1)],2, sd)
  }

  for(i in 1:length(means) ){
    lines(means[[i]], lwd=4, col=cols[i])
    lines(means[[i]]+sdevs[[i]], lwd=1, col=cols[i])
    lines(means[[i]]-sdevs[[i]], lwd=2, col=cols[i])
  }
  }

  legend = c(legend, paste(fly.type[k], c("before", "during", "after") ))
  lty=c(lty, rep(k, 3) )
}

legend(x="topright", col=rep(cols, length(Elist)), legend = legend,
    bty="n", lwd=3, lty = lty) 
