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
WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles, nRuns=1)

#extract segments from features
Slist <- extract.all.segments(WT, k=20)

#do joint procrustes analysis
gpa.joint <- procrustes.analysis.joint(Slist)

#plot first and second principal segment
procrustes.plot(gpa.joint, pc=1, factor=1)
procrustes.plot(gpa.joint, pc=2, factor=3)

#extract the pca transformed coordinate sof the segments for each subject
Z <- procrustes.extract.rawscores.joint(gpa.joint, Slist)

#extracct all sgements based on odor condition
O <- extract.condition.odor.all.list(Z, Slist) 

  
trp <- pairwise.transport(Xin = O, eps=0.05, scale=-1, d=3, store.plan=F)






