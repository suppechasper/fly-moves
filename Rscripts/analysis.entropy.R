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

E1list <- list()
E2list <- list()
E3list <- list()
for(k in 1:length(flies)){
 source(flies[[k]])

  Flist <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles, nRuns=10)

  Slist <- extract.all.segments(Flist, k=1)

  Z <- list()
  for( i in 1:length(Slist) ){
   #Z[[i]] = matrix( Slist[[i]]$lengths, ncol=1) 
   #Z[[i]] = matrix(abs(Slist[[i]]$curvatureMean), ncol=1) 
   Z[[i]] = scale( cbind( Slist[[i]]$lengths, Slist[[i]]$curvatureMean ) )
  }

  O <- list()
  for(i in 1:length(Z) ){
    O[[i]] <- as.matrix( extract.condition.odor2( Z[[i]], Slist[[i]]$C ) )
  }


  E1 <- list()
  E2 <- list()
  E3 <- list()
  for(i in 1:length(O[[1]]) ){
    E1[[i]] <- matrix(0, ncol=scales, nrow=length(O) )
    E2[[i]] <- matrix(0, ncol=scales, nrow=length(O) )
    E3[[i]] <- matrix(0, ncol=scales, nrow=length(O) )
    for(j in 1:length(O) ){
      if(length(O[[j]][[i]]) < 1000){
        E1[[i]][j, ] = NA
        E2[[i]][j, ] = NA
        E3[[i]][j, ] = NA
      }
      else{
        tmp <- multiscale.entropy( O[[j]][[i]], 2, scales, 0.25)
        E1[[i]][j, ] <- tmp[1, ]
        E2[[i]][j, ] <- tmp[2, ]
        E3[[i]][j, ] <- tmp[3, ]
      } 
    }
  }

  E1list[[k]] = E1
  E2list[[k]] = E2
  E3list[[k]] = E3
}


plot(NA, xlim=c(1, scales), ylim= c(0.3, 1.8), bty="n", xlab="scale", ylab="sample entropy" )  
cols = brewer.pal(n=3, "Dark2")
legend=c()
lty=c()

Elist <- E3list
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
