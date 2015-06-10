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
if(F){
WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles, nRuns=10)
}

#extract segments from features
Slist <- extract.all.segments(WT, k=seg.len, off=seg.len/2, rm.jumps=F)

Flist <- extract.all.segment.features(Slist)

library(RColorBrewer)
ramp = colorRamp(c("magenta", "blue"))
cols = rep( c(rgb(t(col2rgb("gold2")/255)), rgb(ramp(seq(0,1, length.out=3) )/255),
    rgb(t(col2rgb("gray"))/255) ), 30)
colsa = paste(cols, "66", sep="") 


X <- c()
C <- c()

for( f in Flist){
  X = rbind(X, f$F)
  C = rbind(C, f$C)
}
X = scale(X)

library(Rtsne)
tsne = Rtsne(X=X, preplexity=25, verbose=T)

symbols = c(1, 2, 3)
plot(tsne$Y, col=cols[C$max.odor], pch=symbols[C$max.pos], xlab="x",
    ylab="y" , asp=1)


library(focus)


#setup and run visualization
data <- focus.create.matrix.data(X)
dataTsne <- focus.create.matrix.data.proxy(data, tsne$Y)

proj2 <- focus.create.orthogonal.projection( dataTsne, fixed=F, hyperbolic=T )
pdel2 = focus.add.projection.display(proj2, 0, 0.5, 0.33, 0.5)

for( r in c(0.25, 0.5, 1.0, 2,0) ){
  focus.add.center.shadow.background(pdel2, r)
}
focus.add.circle.background(pdel2, 1.0)
#focus.add.circle.background(pdel3, 1.0)
#focus.add.circle.background(pdel4, 1.0)


focus.add.profile.display(data, 0.66, 0.0, 0.25, 0.25)

rgbCols = col2rgb(cols)
focus.set.group.colors(data, rgbCols[1, ], rgbCols[2, ], rgbCols[3, ])

focus.set.default.color(data, 0.65, 0, 0)

 
focus.start()

groups = as.factor(C$max.odor)
for(i  in 1:length(levels(groups)) ){
  focus.set.group(data, i-1, which(groups == levels(groups)[i] ) )
}
