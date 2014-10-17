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

flies = c("../Rscripts/WT_ACV0-files.R", "../Rscripts/Orco_ACV0-files.R", "../Rscripts/IR8a1_ACV0-files.R")
fly.type = c("WT", "Orco", "IR8a1")
delay=1
standardize=TRUE

O.all <- list()
W.all <- list()
O.names <- c()


for(i in 1:length(flies)){
 source(flies[[i]])

  #extract features
  WT <- extract.all.features.expected(xyFiles, innerRimFiles, outerRimFiles,
    nRuns=10, std=0.1, lT=5, uT=-1)

  #extract segments from features
   Slist <- extract.all.segments(WT, k=delay)

  Z <- list()
  W <- list()
  for( j in 1:length(Slist) ){
    Z[[j]] = cbind(Slist[[j]]$curvatureMean, Slist[[j]]$lengths)
    if(standardize){
      Z[[j]] = scale( Z[[j]] )
    }
    W[[j]] <- rep(1/length(Z[[j]]), length(Z[[j]]))
  }

  #extracct all sgements based on odor condition
  O <- extract.condition.odor.in.out.all(Z, Slist) 
  names = names( O )
  nperiods = length( names );

  O.all = c(O.all,  O)
  W.all = c(W.all,  W)
  O.names = c( O.names, sprintf("%s-%s", fly.type[i], names) )


}





trp <- pairwise.transport(Xin = O.all, eps=0.0001, scale=-1, d=2, store.plan=T,
    p=2, lambda=0, oType=26, weight=W.all)


#Plot some transport plans
pid = matrix(NA, nrow=length(O.names), ncol=length(O.names) )

for(i in 1:length(trp$plans)){
  pid[ trp$plans[[i]]$i, trp$plans[[i]]$j ] =  i
  pid[ trp$plans[[i]]$j, trp$plans[[i]]$i ] =  i
} 




library(RColorBrewer)
cols = brewer.pal(name="Dark2", n=3)
xs = seq(-3.2, 3.2, length.out=16)
ys= seq(0,5, length.out=16)

for(tp  in c(T, F)){
  for(uc in c(T, F)){

    transparency = tp
      useCost = uc
      lwd= 16
      if(transparency){
        lwd=6
      }



      for(i in 1:(nperiods-1)){
        for(j in (i+1):nperiods){  


            index1= pid[i, j]
            index2= pid[ i + nperiods, j + nperiods   ]
            index3= pid[ i+2*nperiods, j + 2*nperiods ]

            map1 = trp$plans[[ index1 ]]$trp$map[[ length(trp$plans[[index1]]$trp$cost) ]]
            map2 = trp$plans[[ index2 ]]$trp$map[[ length(trp$plans[[index2]]$trp$cost) ]]
            map3 = trp$plans[[ index3 ]]$trp$map[[ length(trp$plans[[index3]]$trp$cost) ]]

            maxW = max( max(map1[,3]), max(map2[,3]), max(map3[,3]) )
            maxC = max( max(map1[,3]*map1[,4]), max(map2[,3]*map2[,4]),
                        max(map3[,3]*map3[,4]) )

            dev.new(width=10, height=10)
            tmp = plot.binnned.transport.map( trp$plans[[index1]]$trp, xs, ys, col=cols[1],
                lwd=lwd, cex.axis=1.5, cex.lab=1.5, useTransparancy=transparency,
                useCost=useCost, maxW=maxW, maxC=maxC )

            tmp = plot.binnned.transport.map( trp$plans[[index2]]$trp, xs, ys, col=cols[2],
                add=T, lwd=lwd, cex.axis=1.5, cex.lab=1.5, useTransparancy=transparency,
                useCost=useCost, maxW=maxW, maxC=maxC )

            tmp = plot.binnned.transport.map( trp$plans[[index3]]$trp, xs, ys, col=cols[3], 
                add=T, lwd=lwd, cex.axis=1.5, cex.lab=1.5, useTransparancy=transparency, 
                useCost=useCost, maxW=maxW, maxC=maxC )

            legend( x="topleft", col=cols, bty="n", legend=fly.type, cex=1.5, lwd=4,
                box.col="gray", bg="white" )

            title( sprintf("%s to %s", names[i], names[j]) )

            if( transparency ){
              if(useCost){
                dev.copy2pdf(file=
                    sprintf("distributions-std-%d-delay-%d-trptc-%s-to-%s.pdf",
                      standardize, delay, names[i], names[j]) )
              }
              else{
                dev.copy2pdf(file= sprintf("distributions-std-%d-delay-%d-trpt-%s-to-%s.pdf", standardize,  delay,
                      names[i], names[j]) )
              } 
            }
            else{
              if(useCost){
                dev.copy2pdf(file= sprintf("distributions-std-%d-delay-%d-trpc-%s-to-%s.pdf", standardize,  delay,
                      names[i], names[j]) ) 
              }
              else{  
                dev.copy2pdf(file= sprintf("distributions-std-%d-delay-%d-trp-%s-to-%s.pdf", standardize,  delay,
                      names[i], names[j]) ) 
              }
            }

            dev.off()
        }
      }

  }
}


#plot distributions

density = list()
library(MASS)
for(i in 1:length(O.all) ){

  density[[i]] = kde2d(x = O.all[[i]][,1], y=O.all[[i]][,2], n=100, lims=c(-3.2, 3.2, 0, 5) )
 
# density[[i]]$z = density[[i]]$z/sum(density[[i]]$z)

            
  dev.new(width=10, height=10)
  
  a = seq(0.01, max(density[[i]]$z)+0.031, by=0.03)
  cols = rgb( colorRamp(brewer.pal(name="Oranges", 9))(seq(0, 1, length.out=length(a))) / 255 )
#contour(density[[i]], levels=a, add=T, col="#00000050", lwd=2)
  filled.contour(density[[i]], levels=a, col=cols, xlab="Curvature", ylab="Speed",
      cex.lab=1.5, bty="n", cex.axis=1.5)
  title( sprintf("Distribution %s", O.names[[i]]))

  dev.copy2pdf(file= sprintf("distributions-std-%d-delay-%d-%s.pdf", standardize,  delay,  O.names[[i]] ) )
  dev.off()

}


#plot differenec distributions
for(i in 1:(length(O.all)-1) ){
  
  d = density[[i]] 
  for(j in (i+1):length(O.all)){
            
    dev.new(width=10, height=10)
    d = density[[i]] 
    d$z = d$z - density[[j]]$z
  
    b = max(abs(min(d$z)), max(d$z))
    a = seq(-b*1.06, b*1.06, length.out=22)
    cols = rgb( colorRamp( rev( brewer.pal("PuOr", n=9)))(seq(0, 1,
            length.out=length(a)-1)) / 255 )
    filled.contour(d, levels=a, col=cols, xlab="Curvature", ylab="Speed",
        cex.lab=1.5, bty="n", cex.axis=1.5)
#contour(d, levels=a, add=T, col="#00000050", lwd=2)
    title( sprintf("Distribution %s minus %s", O.names[[i]], O.names[[j]]))

    dev.copy2pdf(file= sprintf("distributions-std-%d-delay-%d-%s-minus-%s.pdf", standardize,  delay,  O.names[[i]], O.names[[j]] ) )
  
    dev.off()
  } 
      
}





#------------old stuff

if(F){
for(i in 1:length(trp$plans)){

  plan = trp$plans[[i]]$trp
  for(s in 1:length(plan$cost)){
    multiscale.transport.plot.map(plan, s, pointAlpha=0, arrows=T,
        mapAlpha=1, arrow.angle=10, arrow.length=0.1, lwd=2, xlab="curvature",
        ylab="speed")
    title( sprintf("%s to %s scale %i",
          O.names[trp$plans[[i]]$i], O.names[trp$plans[[i]]$j], s ) )
    dev.copy2pdf(file= sprintf("distributions-std-%d-delay-%d-trp-%s-%s-scale-%i.pdf",
          standardize,  delay, O.names[trp$plans[[i]]$i], O.names[trp$plans[[i]]$j], s ) ) 
  }

}
}




