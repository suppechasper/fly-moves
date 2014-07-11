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
source("../Rscripts/WT_ACV0-files.R")


#extract features
WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
  nRuns=10, lT=5, uT=-1)

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
O <- extract.condition.odor.all(Z, Slist) 



procrustes.qq.plot(gpa.joint, pc=3, factor=3, O)


#pairwise t/f tests
TPlist <- list()
Tlist <- list()

FPlist <- list()
Flist <- list()
for(pc in 1:30){
  TPlist[[pc]] = matrix(0, nrow=length(O), ncol=length(O) )
  Tlist[[pc]] = matrix(0, nrow=length(O), ncol=length(O) )
  FPlist[[pc]] = matrix(0, nrow=length(O), ncol=length(O) )
  Flist[[pc]] = matrix(0, nrow=length(O), ncol=length(O) )
}


for(i in 1:(length(O)-1) ){
  for(j in (i+1):length(O) ){
    for(pc in 1:30){
       
       t = t.test(O[[i]][,pc], O[[j]][,pc])
       Tlist[[pc]][i, j] = t$estimate[1] -t$estimate[2]
       Tlist[[pc]][j, i] = t$estimate[1] -t$estimate[2]
       TPlist[[pc]][i, j] = t$p.value
       TPlist[[pc]][j, i] = t$p.value
       
       f = var.test(O[[i]][,pc], O[[j]][,pc])
       Flist[[pc]][i, j] = f$estimate
       Flist[[pc]][j, i] = f$estimate
       FPlist[[pc]][i, j] = f$p.value
       FPlist[[pc]][j, i] = f$p.value

    }
  }
}

