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
  nRuns=10, std=0.1, lT=5, uT=0.25)

#extract segments from features
Slist <- extract.all.segments(WT, k=1)

Z <- list()
for( i in 1:length(Slist) ){
  Z[[i]] = cbind(Slist[[i]]$curvatureMean, Slist[[i]]$lengths)
}

#extracct all sgements based on odor condition
O <- extract.condition.odor.all(Z, Slist) 

  
trp <- pairwise.transport(Xin = O, eps=0.0001, scale=-1, d=2, store.plan=T, p=1)

multiscale.transport.plot.map(trp$plans[[2]]$trp, 8, pointAlpha=0, arrows=T,
    mapAlpha=0.7, arrow.angle=10, arrow.length=0.1, lwd=2)




