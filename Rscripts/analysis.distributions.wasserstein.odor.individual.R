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
WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles, nRuns=10)

#extract segments from features
Slist <- extract.all.segments(WT, k=1)

Z = list()
for( i in 1:length(Slist) ){
  Z[[i]] = cbind(Slist[[i]]$curvatureMean, Slist[[i]]$lengths)
}

D <- c()
O.trp <- list()
paths <- list()

for(i in 1:length(Z)){
  #extracct all sgements based on odor condition
  O <- extract.condition.odor(Z[[i]], Slist[[i]]$C) 

  #put all oders in a list
  Olist <- list(O$Xbefore2, O$Xbefore2, O$Xduring1, O$Xduring2, O$Xafter1,
    O$Xafter2)


  trp <- pairwise.transport(Xin = Olist, eps=0.0001, scale=-1, store.plan=T)

  P <- list()
  for(j in 1:length(trp$plans)){
    P[[j]] <- transport.extract.paths(trp$plans[[j]]$trp,
        length(trp$plans[[j]]$trp$cost) )
  }
  paths[[i]] = P
  O.trp[[i]] = trp;
  
  D = rbind( scale(as.vector(trp$D)) )
}








