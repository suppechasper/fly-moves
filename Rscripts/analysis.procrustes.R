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
seg.len = 20
flies = c("../Rscripts/WT_ACV0-files.R", "../Rscripts/Orco_ACV0-files.R", "../Rscripts/IR8a1_ACV0-files.R")
fly.type = c("WT", "Orco", "IR8a1")

for(k in 1:length(flies)){
 source(flies[[k]])

#extract features
WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
  nRuns=10, lT=5, uT=0.1)

#extract segments from features
Slist <- extract.all.segments(WT, k=seg.len)

#do joint procrustes analysis
gpa.joint <- procrustes.analysis.joint(Slist)

#plot first and second principal segment
#procrustes.plot(gpa.joint, pc=1, factor=1)
#procrustes.plot(gpa.joint, pc=2, factor=3)

#extract the pca transformed coordinate sof the segments for each subject
Z <- procrustes.extract.rawscores.joint(gpa.joint, Slist, 1:10)
  
#extracct all sgements based on odor condition
O <- extract.condition.odor.all(Z, Slist) 
O.names = names(O)


for(i in 1:10){
  f = 3;
  if(i==1){
    f=1
  }
  procrustes.qq.plot(gpa.joint, pc=i, factor=f, O, 10000)
  dev.copy2pdf( file=sprintf("segments-%s-len-%d-qqplot-pc-%d.pdf", fly.type[k], seg.len, i) )
}











if(F){


trps <- pairwise.transport(Xin = O, eps=0.01, scale=-1, rFactor=0.5, d=3,
    store.plan=T, p=2, lambda=0.05, oType=26, split=2, sType=0)

  
for(i in 1:length(trps$plans)){

  plan = trps$plans[[i]]$trp

  for(s in 2:length(plan$cost)){

    par(mar=c(5,5,5,5))
    layout( matrix(c(rep(1, 9), 2:7), nrow=3) )

    pc = plot.pc.transport.map(plan, s)

    title( sprintf("%s to %s scale %i",
          O.names[trps$plans[[i]]$i], O.names[trps$plans[[i]]$j], s ) )

    d = dim( gpa.joint$mshape )

    dir = gpa.joint$pcar[,1:10] %*% pc$rotation[,1]
    dir = dir/sqrt( sum(dir^2) )
    dir = 3*pc$sdev[1]*matrix(dir, nrow=d[1])
    pc.segment.plot(gpa.joint$mshape, dir, 3, 1)
    
    dir = gpa.joint$pcar[,1:10] %*% pc$rotation[,2]
    dir = dir/sqrt( sum(dir^2) )
    dir = 3*pc$sdev[2]*matrix(dir, nrow=d[1])
    pc.segment.plot(gpa.joint$mshape, dir, 3, 2)

    dev.copy2pdf(file= sprintf("segments-%s-len-%d-trp-%s-%s-scale-%i.pdf",
          fly.type[k], seg.len, O.names[trps$plans[[i]]$i], O.names[trps$plans[[i]]$j],
          s )) 
 
  }

}

}
 










if(F){

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
}



}
