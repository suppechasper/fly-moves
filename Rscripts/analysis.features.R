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

#extract segments from features
Slist <- extract.all.segments(WT, k=seg.len, off=seg.len/2)

Flist <- extract.all.segment.features(Slist)

library(RColorBrewer)
ramp = colorRamp(c("magenta", "blue"))
cols = rep( c(rgb(t(col2rgb("gold2")/255)), rgb(ramp(seq(0,1, length.out=3) )/255),
    rgb(t(col2rgb("black"))/255) ), 30)
colsa = paste(cols, "66", sep="") 
