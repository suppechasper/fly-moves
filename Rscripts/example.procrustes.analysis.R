source("extract.features.R")
source("procrustes.R")
source("delay.reconstruction.R")

#to install packages uncomment:
#install.packages("shapes")
#install.packages("abind")
#install.packages("circular")

#set path to data files
#list of file names to be read is set in extract.features.R
path = "../../WT_ACV0/orig-data/"


#extract features
F <- extract.all.features.expected()

#extract segments from features
Slist <- extract.all.segments(F, k=20)

#do joint procrustes analysis
gpa.joint <- procrustes.analysis.joint(Slist)

#plot first and second principal segment
procrustes.plot(gpa.joint, pc=1, factor=1)
procrustes.plot(gpa.joint, pc=2, factor=3)

#extract the pca transformed coordinate sof the segments for each subject
Z <- procrustes.extract.rawscores.joint(gpa.joint, Slist)
  
#extracct all sgements based on odor condition
O <- extract.condition.odor.all(Z, Slist) 


#student t test of difference in pc
#not particularly helpful...
pc = 1
t.test(O$Xbefore1[,pc], O$Xbefore2[,pc])
t.test(O$Xbefore1[,pc], O$Xduring1[,pc])
t.test(O$Xbefore1[,pc], O$Xduring2[,pc])
t.test(O$Xbefore1[,pc], O$Xafter1[,pc])
t.test(O$Xbefore1[,pc], O$Xafter2[,pc])

t.test(O$Xbefore2[,pc], O$Xduring1[,pc])
t.test(O$Xbefore2[,pc], O$Xduring2[,pc])
t.test(O$Xbefore2[,pc], O$Xafter1[,pc])
t.test(O$Xbefore2[,pc], O$Xafter2[,pc])

t.test(O$Xduring1[,pc], O$Xduring2[,pc])
t.test(O$Xduring1[,pc], O$Xafter1[,pc])
t.test(O$Xduring1[,pc], O$Xafter2[,pc])

t.test(O$Xduring2[,pc], O$Xafter1[,pc])
t.test(O$Xduring2[,pc], O$Xafter2[,pc])

t.test(O$Xafter1[,pc], O$Xafter2[,pc])


