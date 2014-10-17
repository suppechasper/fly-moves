source("../Rscripts/extract.features.R")
source("../Rscripts/procrustes.R")
source("../Rscripts/delay.reconstruction.R")

#to install packages uncomment:
#install.packages("shapes")
#install.packages("abind")
#install.packages("circular")

#load WT file names, adjust path in WT_ACV0-files.R tp point to the correct
#directoty. To convert .mat files (rim information files) to csv use
#rim-mat-2-csv.m in OctaveScripts (should be compatible with Matlab

flies = c("../Rscripts/WT_ACV0-files.R", "../Rscripts/Orco_ACV0-files.R", "../Rscripts/IR8a1_ACV0-files.R")
fly.type = c("WT", "Orco", "IR8a1")
delay=20
standardize=TRUE

O <- list()
O.names <- list()


for(i in 1:length(flies)){
 source(flies[[i]])

  #extract features
  WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
    nRuns=10, lT=5, uT=-1)

  #extract segments from features
   Slist <- extract.all.segments(WT, k=delay)

  Z <- list()
  W <- list()
  for( j in 1:length(Slist) ){
    Z[[j]] = cbind(  
                   rowMeans( Slist[[j]]$curvatureMean ), 
                   sqrt( apply(Slist[[j]]$curvatureMean, 1, var) ), 
                   rowMeans(Slist[[j]]$lengths), 
                   sqrt( apply(Slist[[j]]$lengths, 1, var) )
                 )
    
    if(standardize){
      Z[[j]] = scale( Z[[j]] )
    }
    W[[j]] <- rep(1/length(Z[[j]]), length(Z[[j]]))
  }

  #extracct all sgements based on odor condition
  O[[i]] <- extract.condition.odor.in.out.all(Z, Slist) 
  names = names( O[[i]] )
  nperiods = length( names );

  O.names[[i]] = sprintf("%s-%s", fly.type[i], names) 


}


#plot per features densities
library(RColorBrewer)
dev.new(width=10, height=10)

colsTmp = brewer.pal(name="Set1", n=length(O))
cols = c()
lty = c()
for(i in 1:length(colsTmp) ){
  cols =c(cols, rep(colsTmp[i], 3) )
  lty = c(lty, 1:3) 
} 
colsa <- sprintf("%s%s", cols, "55")
for(i in 1:4){
  f=3
  if(i==1){
    f=1
  }
  for(k in 1:length(O) ){
    procrustes.density.plot2(O[[k]], pc=i, factor=f, col=cols, lty=lty,
        legend=names(O[[k]]), lwd=2 )
    dev.copy2pdf( file=sprintf("segments-len-%.3d-fly-%s-density-plot-pc-%.3d.pdf",
          delay, fly.type[k], i)  )
  }
}




