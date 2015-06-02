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


#extract features
WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles, nRuns=10)

#extract segments from features
Slist <- extract.all.segments(WT, k=1)

Z = list()
for( i in 1:length(Slist) ){
  Z[[i]] =  cbind(Slist[[i]]$curvatureMean, Slist[[i]]$lengths, Slist[[i]]$orientation )  
}

D <- c()
Olist <- list()
paths <- list()


for(i in 1:length(Z)){
  #extract all sgements based on odor condition
  O <- extract.condition.odor(Z[[i]], Slist[[i]]$C) 
  #put all odors in a list
  Olist <- c(Olist, O)
}


Zbefore = c()
for(i in seq(1, length(Olist), by=5 )){
  Zbefore = rbind( Zbefore, Olist[[i]] )
}
aM = colMeans(Zbefore )
aV = diag( var(Zbefore) )



saOlist = list()
for(i in 1:length(Olist) ){
    saOlist[[i]] = Olist[[i]]
    for(k in 1:ncol(saOlist[[i]]) ){
        saOlist[[i]][, k] = saOlist[[i]][, k] - aM[k]
        saOlist[[i]][, k] = saOlist[[i]][, k] / sqrt(aV[k])
    }
}



V = matrix( 0, nrow=length(saOlist), ncol=3 )
M = matrix( 0, nrow=length(saOlist), ncol=3 )
for(i in 1:length(saOlist) ){
  M[i, ] = colMeans(saOlist[[i]] )
  V[i, ] = diag( var(saOlist[[i]]) )
}



sOlist = list()
for(i in seq(1, length(saOlist), by=5 )){
  for(j in 0:4){
    sOlist[[i+j]] = saOlist[[i+j]]
    for(k in 1:ncol(sOlist[[i+j]]) ){
        sOlist[[i+j]][, k] = sOlist[[i+j]][, k] - M[i,k]
        sOlist[[i+j]][, k] = sOlist[[i+j]][, k] / sqrt(V[i,k])
    }
  }
}

Vs = matrix( 0, nrow=length(sOlist), ncol=3 )
Ms = matrix( 0, nrow=length(sOlist), ncol=3 )
for(i in 1:length(sOlist) ){
  Ms[i, ] = colMeans(sOlist[[i]] )
  Vs[i, ] = diag( var(sOlist[[i]]) )
}

#color by odor
ramp = colorRamp(c("magenta", "blue"))
colsAll = rep( c("gold2", rgb(ramp(seq(0,1, length.out=3) )/255), "black"), 30)
cols = colsAll


plot(M[, 1:2], col = cols, bty="n", xlab="mean( curvature )", ylab= "mean( speed )")
dev.copy2eps(file = "mean-12.eps")
plot(M[, c(1,3)], col = cols, bty="n", xlab="mean( curvature )", ylab= "mean( orientation )")
dev.copy2eps(file = "mean-13.eps")
plot(V[, 1:2], col = cols, bty="n", xlab="var( curvature )", ylab= "var( speed )")
dev.copy2eps(file = "var-12.eps")
plot(V[, c(1,3)], col = cols, bty="n", xlab="var( curvature )", ylab= "var( orientation )")
dev.copy2eps(file = "var-13.eps")

plot(Ms[, 1:2], col = cols, bty="n", xlab="mean( curvature )", ylab= "mean( speed )")
dev.copy2eps(file = "mean-scaled-12.eps")
plot(Ms[, c(1,3)], col = cols, bty="n", xlab="mean( curvature )", ylab= "mean( orientation )")
dev.copy2eps(file = "mean-scaled-13.eps")
plot(Vs[, 1:2], col = cols, bty="n", xlab="var( curvature )", ylab= "var( speed )")
dev.copy2eps(file = "var-scaled-12.eps")
plot(Vs[, c(1,3)], col = cols, bty="n", xlab="var( curvature )", ylab= "var( orientation )")
dev.copy2eps(file = "var-scaled-13.eps")

