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

prefix="ACV0"
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
O.names = names(O)
  
trp <- pairwise.transport(Xin = O, eps=0.0001, scale=-1, d=2, store.plan=T, p=1)

for(i in 1:length(trp$plans)){
  
  plan = trp$plans[[i]]$trp
  for(s in 1:length(plan$cost)){
    multiscale.transport.plot.map(plans, s, pointAlpha=0, arrows=T,
        mapAlpha=1, arrow.angle=10, arrow.length=0.1, lwd=2, xlab="curvature",
        ylab="speed")
    title( sprintf("%s to %s scale %i",
          O.names[trp$plans[[i]]$i], O.names[trp$plans[[i]]$j], s ) )
  }

}




