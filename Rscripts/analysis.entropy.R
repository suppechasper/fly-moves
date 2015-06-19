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
seg.len = 20

#extract features
WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles, nRuns=10)

Slist <- extract.all.segments(WT, k=1)

Z <- list()
for( i in 1:length(Slist) ){
  Z[[i]] = matrix(Slist[[i]]$lengths, ncol=1) 
}

O <- list()
for(i in 1:length(Z) ){
  O[[i]] <- extract.condition.odor2( Z[[i]], Slist[[i]]$C ) 
}


E <- list()
for(i in 1:length(O) ){
  E[[i]] <- matrix(0, ncol=20, nrow=length(O[[i]]) )
  for(j in 1:length(O[[i]]) ){
    E[[i]][j, ] <- multiscale.entropy(matrix(O[[i]][[j]], ncol=1), 2, 20, 0.15)[1, ] 
  }
}
