source("../Rscripts/extract.features.R")
source("../Rscripts/procrustes.R")

#to install packages uncomment:
#install.packages("shapes")
#install.packages("abind")
#install.packages("circular")
#install.packages("KernSMooth")

#load WT file names, adjust path in WT_ACV0-files.R tp point to the correct
#directoty. To convert .mat files (rim information files) to csv use
#rim-mat-2-csv.m in OctaveScripts (should be compatible with Matlab


seg.lens = c(5, 10, 20, 40, 80, 160)
lT = 5
uT = -1
for(ii in 1:length(seg.lens) ){

seg.len = seg.lens[ii]
offset <- floor(seg.len/2)

flies <- c("../Rscripts/WT_ACV0-files.R", "../Rscripts/Orco_ACV0-files.R", "../Rscripts/IR8a1_ACV0-files.R")
fly.type <- c("WT", "Orco", "IR8a1")

Slist <- c()
fly.type.id <- c()
fly.ids <- c()
for(k in 1:length(flies)){

  source(flies[[k]])

  #extract features
  WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
    nRuns=10, lT=lT, uT=uT)

  #extract segments from features
  Stmp <- extract.all.segments(WT, k=seg.len, offset=offset)
 
  fly.ids <- c(fly.ids, rep(k, length(Stmp)))

  for(i in 1:length(Stmp)){
    fly.type.id <- c(fly.type.id, rep(k, nrow(Stmp[[i]]$C)))
  }
  Slist <- c(Slist, Stmp)
}

times <- c()
for(i in 1:length(Slist) ){
  times <- c(times, Slist[[i]]$C$timeToOdor)
}



#do joint procrustes analysis
gpa.joint <- procrustes.analysis.joint(Slist)
save(gpa.joint,
    file=sprintf("gpa-joint-seegments-len-%.3d-lt-%.5f-ut-%.5f.Rdata", seg.len,
      lT, uT) )
npcs = min(ncol(gpa.joint$pcar), 20)
#extract the pca transformed coordinates of the segments for each subject
Z <- procrustes.extract.rawscores.joint(gpa.joint, Slist, 1:npcs)
  
nperiods = 0
O <- list()
for(k in 1:3){
    ind = which(fly.ids==k)
    O[[k]] <- extract.condition.odor.in.out.all(Z[ind], Slist[ind]) 
    nperiods = length(O[[k]])
}


#par(mar=c(5,5,5,5))

#plot per pcs spectrums

dev.new(width=40, height=10)
for(i in 1:npcs ){
  procrustes.plot.pcs(gpa.joint, i,  O) 
  dev.copy2pdf( file=sprintf("segments-len-%.3d-spectrum-pc-%.3d.pdf", seg.len, i) )
}
dev.off()
  



#plot spectrums
dev.new(width=20, height=10)
Vtotal = c()
Vfly= c()
Vodor = c()
for(pc in 1:10){
    se = c()
    counts = c()
    nperiods = 0
    for(i in 1:length(O)){
      nperiods = length(O[[i]])
        for(k in 1:length(O[[i]]) ){
          se = c(se, sum(O[[i]][[k]][,pc]^2))
          counts = c(counts, nrow(O[[i]][[k]]) )
        }
    } 
    Vtotal = c(Vtotal, sum(se)/(sum(counts)-1))
    
    Vperfly = c()
    for(i in seq(1, length(se), by=nperiods)){
      Vperfly = c(Vperfly,
      sum(se[i:(i+nperiods-1)])/(sum(counts[i:(i+nperiods-1)])-1) )
      
    }
    Vfly = rbind(Vfly, Vperfly) 

    Vperflyperodor = se/(counts-1)
    Vodor = rbind(Vodor, Vperflyperodor) 
}

colnames(Vfly)  = NULL
rownames(Vfly) = NULL
colnames(Vodor) = NULL
rownames(Vodor) = NULL

barplot(Vtotal)
dev.copy2pdf( file=sprintf("segments-len-%.3d-spectrum-total.pdf", seg.len) )

library(RColorBrewer)
colsTmp = brewer.pal(name="Dark2", n=length(O))


for(k in 1:ncol(Vfly)){
  barplot(Vfly[, k], col=colsTmp[k] )
  dev.copy2pdf( file=sprintf("segments-len-%.3d-spectrum-fly-%s.pdf", seg.len,
        fly.type[k] ) )
  for(i in 1:nperiods ){
    barplot(Vodor[,i+((k-1)*nperiods)], col=colsTmp[k] )
    dev.copy2pdf( file=sprintf("segments-len-%.3d-spectrum-fly-%s-odor-%s.pdf",
          seg.len, fly.type[k], names(O[[k]])[i] ) ) 
  }


}



#plot per pc densities
dev.off()
dev.new(width=10, height=10)

colsTmp = brewer.pal(name="Set1", n=length(O))
cols = c()
lty = c()
for(i in 1:length(colsTmp) ){
  cols =c(cols, rep(colsTmp[i], 3) )
  lty = c(lty, 1:3) 
} 
colsa <- sprintf("%s%s", cols, "55")
for(i in 1:npcs){
  f=3
  if(i==1){
    f=1
  }
  for(k in 1:length(O) ){
    procrustes.density.plot2(O[[k]], pc=i, factor=f, col=cols, lty=lty, legend=names(O[[k]]), lwd=2 )
    dev.copy2pdf( file=sprintf("segments-len-%.3d-fly-%s-density-plot-pc-%.3d.pdf",
          seg.len, fly.type[k], i)  )
  }
}








}





#-----------older stuff
if(F){
mean = gpa.joint$mshape
d = dim( gpa.joint$mshape )
par(mar=c(5,5,5,5))
pca <- list()
for(k in 1:3){

    pca[[k]] = list()
    for( i in 1:length(O[[k]])){
       pca[[k]][[i]] = prcomp(O[[k]][[i]])
    }
}

for(pc in 1:10){
      factor = 3;
  if(i==1){
    factor=1
  }

    layout( matrix(1:18, nrow=3) )


    for(i in 1:length(pca[[1]])){

      for(s in c(-1,0,1)){
      for(k in 1:length(pca)){

        dir = gpa.joint$pcar[,1:10] %*% pca[[k]][[i]]$rotation[,pc]
        dir = dir/sqrt( sum(dir^2) )
        dir = 3*pca[[k]][[i]]$sdev[pc]*matrix(dir, nrow=d[1])  
        n = nrow(mean)
        xlim = range(mean[,1]) + range(dir[,1])
        ylim = range(mean[,2]) + range(dir[,2])
        if(k==1){  
          plot( mean +s*dir, pch=19, xlim=xlim, ylim=ylim, asp=1, bty="n", xlab="",
                ylab="", cex.axis=2, cex=2, col=colsa[k], type="l", lwd=5)
          points((mean +s* dir)[2:n, ], pch=19, col=colsa[k],cex=2)
          title(sprintf("Mean %d sd * PC %d", s*factor, pc), cex.main=2)
        }
        else{
          lines( mean + s* dir, pch=19, xlim=xlim, ylim=ylim, asp=1, bty="n", xlab="",
                 ylab="", cex.axis=2, cex=2, col=colsa[k], type="l", lwd=5)
          points((mean +s* dir)[2:n, ], pch=19, col=colsa[k],cex=2)
        }
      }
      
      }

    }


    dev.copy2pdf(file= sprintf("segments-len-%d-pc-%d.pdf", seg.len, pc)) 

}



}














##lda test
if(F){
X = rbind(O[[3]][[1]], O[[2]][[1]])
g = c( rep("before", nrow(O[[3]][[1]])), rep("during", nrow(O[[2]][[1]])))

q <- lda(x=as.matrix(X), as.factor(g))

p <- q$scaling / sqrt(sum(q$scaling^2))
pq <- X %*% p
spq = sd(pq)

d = dim( gpa.joint$mshape )
    
dir = gpa.joint$pcar[,1:10] %*% q$scaling
dir = dir/sqrt( sum(dir^2) )
dir = spq*matrix(dir, nrow=d[1])  
layout( t(matrix(c(1:3, 4,4,4), nrow=3)) )
pc.segment.plot(gpa.joint$mshape, dir, 1, 1)
d1 <- bkde(x=pq[g=="before"])
d2 <- bkde(x=pq[g=="during"])

plot(d1, type="l", bty="n", lwd=3)
lines(d2, lwd=3, col="red")
}











####
if(F){

library(RColorBrewer)
cols = brewer.pal(name="Dark2", n=3)

for(i in 1:10){
  f = 3;
  if(i==1){
    f=1
  }

  par(mar=c(5,5,5,5))
  d = dim( gpa.joint$mshape )
  dir = f* gpa.joint$pcasd[i] *matrix(gpa.joint$pcar[,i] , nrow=d[1])
  
  layout( matrix(c(1:3, rep(4,3*3)), nrow=3) )
  pc.segment.plot(gpa.joint$mshape, dir, f, i)


  procrustes.time.plot(gpa.joint, pc=i, factor=f, Z=gpa.joint$rawscores,
      times=times, ncuts=250, bw=250, col="#00000050")
  abline(h=0, col="#00000020", lwd=1)
  abline(v=0, col="#00000020", lwd=1)
  for(k in 1:3){
    ind = which(fly.type.id==k)
    procrustes.time.plot(gpa.joint, pc=i, factor=f, Z=gpa.joint$rawscores[ind, ],
      times=times[ind], ncuts=250, bw=250, col=sprintf("%s50", cols[k]), add=T)
  }

  dev.copy2pdf( file=sprintf("segments-%s-len-%d-time-plot-pc-%d.pdf", fly.type[k], seg.len, i) )
 
  procrustes.qq.plot(gpa.joint, pc=i, factor=f, O, 10000)
}

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
