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
scales = 20

#extract features
WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles, nRuns=10)

Slist <- extract.all.segments(WT, k=1)

Z <- list()
for( i in 1:length(Slist) ){
  Z[[i]] = matrix(abs(Slist[[i]]$curvatureMean), ncol=1) 
}

O <- list()
for(i in 1:length(Z) ){
  O[[i]] <- extract.condition.odor2( Z[[i]], Slist[[i]]$C ) 
}


library(mmsentropy)

E <- list()
for(i in 1:length(O[[1]]) ){
  E[[i]] <- matrix(0, ncol=scales, nrow=length(O) )
  for(j in 1:length(O) ){
    E[[i]][j, ] <- multiscale.entropy(matrix(O[[j]][[i]], ncol=1), 2, scales,
        0.25)[1, ] 
  }
}


means <- list()
sdevs <- list()
for(i in 1:length(E) ){
  means[[i]] = colMeans(E[[i]])
  sdevs[[i]] = apply(E[[i]],2, sd)
}


library(RColorBrewer)

plot(NA, xlim=c(1, scales), ylim= c(min(means[[i]]) - 0.5, max(means[[i]]) + 0.5), bty="n", xlab="scale", ylab="sample entropy" )  
cols = brewer.pal(n=length(means), "Dark2")
for(i in 1:length(means) ){
  lines(means[[i]], lwd=4, col=cols[i])
  lines(means[[i]]+sdevs[[i]], lwd=1, col=cols[i])
  lines(means[[i]]-sdevs[[i]], lwd=2, col=cols[i])
}
legend(x="bottomright", col=cols, legend = c("before", "during", "after"),
    bty="n", lwd=3) 

means <- list()
sdevs <- list()
for(i in 1:length(E) ){
  means[[i]] = colMeans(E[[i]][, 2:scales] - E[[i]][, 1:(scales-1)])
  sdevs[[i]] = apply(E[[i]][, 2:scales] - E[[i]][, 1:(scales-1)],2, sd)
}

dev.new()
plot(NA, xlim=c(1, scales), ylim= c(min(means[[i]]) - 0.5, max(means[[i]]) + 0.5), bty="n", xlab="scale", ylab="sample entropy" )  
for(i in 1:length(means) ){
  lines(means[[i]], lwd=4, col=cols[i])
  lines(means[[i]]+sdevs[[i]], lwd=1, col=cols[i])
  lines(means[[i]]-sdevs[[i]], lwd=2, col=cols[i])
}


