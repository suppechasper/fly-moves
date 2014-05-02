##
#frame odor on:   5400
#frame odor off: 10801
#fly dimensions 25 x 12 pixels
#instantaneous velocity < 3 pixels /s , values abouve 4 p/s probably errorenous
#peak fluctation is 0.15 pixels, sd 0.09 pixels

conditional.samples <- function(data, conds, n=30){

  odor = conds$time > 5400 & conds$time < 10801
  inside = conds$dirim > -25
  atRim = abs(conds$dorim) < 25
  id = conds$id 


  #build  a set of conditional samples
  X <- list()
  
  C <- matrix(0, nrow= n * 18 + 18, ncol=6)
  colnames(C) <- c("odorOn", "odorOff", "innerRim", "middle", "outerRim", "id") 


  index = 1;
  #conditional averages over all flies
  X[[index]] = data
  C[index, ] = 1
  index = index + 1
  
  #
  X[[index]] = data[inside == 1, ]
  C[index, ] = c(1, 1, 1, 0, 0, -1)
  index = index +1
  
  X[[index]] = data[inside == 0 & atRim == 0, ]
  C[index, ] = c(1, 1, 0, 1, 0, -1)
  index = index +1

  X[[index]] = data[ inside == 0, ]
  C[index, ] = c(1, 1, 0, 1, 1, -1)
  index = index +1
  
  X[[index]] = data[ atRim==1, ]
  C[index, ] = c(1, 1, 0, 0, 1, -1)
  index = index + 1

  X[[index]] = data[ atRim==0, ]
  C[index, ] = c(1, 1, 1, 1, 0, -1)
  index = index + 1


  #
  X[[index]] = data[odor == 1 & inside == 1, ]
  C[index, ] = c(1, 0, 1, 0, 0, -1)
  index = index +1
  
  X[[index]] = data[odor == 1 & inside == 0 & atRim == 0, ]
  C[index, ] = c(1, 0, 0, 1, 0, -1)
  index = index +1

  X[[index]] = data[odor == 1 & inside == 0, ]
  C[index, ] = c(1, 0, 0, 1, 1, -1)
  index = index +1
  
  X[[index]] = data[odor==1 & atRim==1, ]
  C[index, ] = c(1, 0, 0, 0, 1, -1)
  index = index +1

  X[[index]] = data[odor==1, ]
  C[index, ] = c(1, 0, 1, 1, 1, -1)
  index = index +1

  X[[index]] = data[ odor==1 & atRim==0, ]
  C[index, ] = c(1, 0, 1, 1, 0, -1)
  index = index + 1


  #
  X[[index]] = data[odor == 0 & inside == 1, ]
  C[index, ] = c(0, 1, 1, 0, 0, -1)
  index = index +1
  
  X[[index]] = data[odor == 0 & inside == 0 & atRim == 0, ]
  C[index, ] = c(0, 1, 0, 1, 0, -1)
  index = index +1

  X[[index]] = data[odor == 0 & inside == 0, ]
  C[index, ] = c(0, 1, 0, 1, 1, -1)
  index = index +1
  
  X[[index]] = data[odor==0 & atRim==1, ]
  C[index, ] = c(0, 1, 0, 0, 1, -1)
  index = index +1

  X[[index]] = data[odor==0, ]
  C[index, ] = c(0, 1, 1, 1, 1, -1)
  index = index +1

  X[[index]] = data[ odor==0 & atRim==0, ]
  C[index, ] = c(0, 1, 1, 1, 0, -1)
  index = index + 1






  #per fly conditions
  for(i in 1:max(conds$id) ){
    print(i)

  X[[index]] = data[id == i, ]
  C[index, ] = c(1, 1, 1, 1, 1, i)
  index = index + 1
  
  #
  X[[index]] = data[inside == 1 & id == i, ]
  C[index, ] = c(1, 1, 1, 0, 0, i)
  index = index +1
  
  X[[index]] = data[inside == 0 & atRim == 0 & id == i, ]
  C[index, ] = c(1, 1, 0, 1, 0, i)
  index = index +1

  X[[index]] = data[ inside == 0 & id == i, ]
  C[index, ] = c(1, 1, 0, 1, 1, i)
  index = index +1
  
  X[[index]] = data[atRim==1 & id == i, ]
  C[index, ] = c(1, 1, 0, 0, 1, i)
  index = index +1
  
  X[[index]] = data[ atRim==0 & id==i, ]
  C[index, ] = c(1, 1, 1, 1, 0, i)
  index = index + 1


  #
  X[[index]] = data[odor == 1 & inside == 1 & id == i, ]
  C[index, ] = c(1, 0, 1, 0, 0, i)
  index = index +1
  
  X[[index]] = data[odor == 1 & inside == 0 & atRim == 0 & id == i, ]
  C[index, ] = c(1, 0, 0, 1, 0, i)
  index = index +1

  X[[index]] = data[odor == 1 & inside == 0 & id == i, ]
  C[index, ] = c(1, 0, 0, 1, 1, i)
  index = index +1
  
  X[[index]] = data[odor==1 & atRim==1 & id == i, ]
  C[index, ] = c(1, 0, 0, 0, 1, i)
  index = index +1

  X[[index]] = data[odor==1 & id == i, ]
  C[index, ] = c(1, 0, 1, 1, 1, i)
  index = index +1

  X[[index]] = data[ odor==1 & atRim==0 & id==i, ]
  C[index, ] = c(1, 0, 1, 1, 0, i)
  index = index + 1


  #
  X[[index]] = data[odor == 0 & inside == 1 & id == i, ]
  C[index, ] = c(0, 1, 1, 0, 0, i)
  index = index +1
  
  X[[index]] = data[odor == 0 & inside == 0 & atRim == 0 & id == i, ]
  C[index, ] = c(0, 1, 0, 1, 0, i)
  index = index +1

  X[[index]] = data[odor == 0 & inside == 0 & id == i, ]
  C[index, ] = c(0, 1, 0, 1, 1, i)
  index = index +1
  
  X[[index]] = data[odor==0 & atRim==1 & id == i, ]
  C[index, ] = c(0, 1, 0, 0, 1, i)
  index = index +1

  X[[index]] = data[odor==0 & id == i, ]
  C[index, ] = c(0, 1, 1, 1, 1, i)
  index = index +1

  X[[index]] = data[ odor==0 & atRim==0 & id==i, ]
  C[index, ] = c(0, 1, 1, 1, 0, i)
  index = index + 1


  }

  n <- unlist( lapply(X, nrow) )

  res = list( samples = X, conditions=C, nSamples = n)

}





conditional.procrustes <- function(S){
  library(shapes)

  C = as.data.frame( S$conditions )

  ind = which(C$odorOn == 1 & C$odorOff == 0 & C$innerRim==1 & C$middle==0 &
      C$outerRim ==0 & C$id == -1 & S$nSamples > 0)
  X <- extract.procrustes.segments.block(S$samples[ind])
  gpa <- procGPA(X, scale=F, reflect=T, proc.output=T, distances=F, pcaoutput=F)
  save(gpa, ind, file="segs-20-gpa-odor-inside.Rdata")

  ind = which( C$middle==1 & C$innerRim==0 & C$outerRim==0 & C$id==-1, S$nSamples > 0)
  X <- extract.procrustes.segments.block(S$samples[ind])
  gpa <- procGPA(X, scale=F, reflect=T, proc.output=T, distances=F, pcaoutput=F)
  save(gpa, ind, file="segs-20-gpa-middle.Rdata")

  ind = which( C$middle==0 & C$innerRim==0 & C$outerRim==1 & C$id==-1 & S$nSamples > 0)
  X <- extract.procrustes.segments.block(S$samples[ind])
  gpa <- procGPA(X, scale=F, reflect=T, proc.output=T, distances=F, pcaoutput=F)
  save(gpa, ind, file="segs-20-gpa-atOuterRim.Rdata")

}





extract.procrustes.segments.block <- function(data){
  library(shapes)
  library(abind)
  
  l = ncol(data[[1]])
  s = l/2
  segs <- c()
  for(seg in data){
     x = abind::abind(seg[,1:s], along=1)
     y = abind::abind(seg[,(s+1):l], along=1)
     z = aperm( abind::abind(x, y, along=3), c(2,3,1) )
     segs = abind::abind(segs, z, along=3)
  }
  segs
}





extract.procrustes.segments <- function(data){
  library(shapes)
  library(abind)
  k <- length( data[[1]]$segsX1 )
  segs <- c()
  for(seg in data){
     x = abind::abind(seg$segsX1, along=1)
     y = abind::abind(seg$segsX2, along=1)
     z = aperm( abind::abind(x, y, along=3), c(2,3,1) ) 
     segs = abind::abind(segs, z, along=3)
  }
  segs
}





pairwise.transport <- function(Xin, eps, scale=-1){
  library(mop)


  gmra <- c()
  indices = c()
  for(i in 1:length(Xin)){
    nr = nrow(Xin[[i]])
    if(!is.null(nr)){
      if( nr > 0){

        gmra <- c(gmra, multiscale.transport.create.ipca(X=Xin[[i]], d=3,
            eps=eps, t=0.9) )
        indices = c(indices, i)
      }
    }
  }
  
  dist = matrix( 0, nrow=length(gmra), ncol=length(gmra) )

  for(i in 1:(length(gmra)-1)){
    for(j in (i+1):length(gmra)){
      trp = multiscale.transport.id(gmra1=gmra[i], gmra2 = gmra[j],
          p=2, rFactor=1, sType=0, scale1=scale, scale2=scale, stpPct=-1,
          oType=26, propFactor=0,  nRefinementIterations = 1 )
     
      dist[i, j] = trp$cost[length(trp$cost)]
      dist[j, i] = dist[i, j]
#save(trp, file=sprintf("%s-%d-%d.Rdata", prefix, i, j))
    }
  }

  list(D = dist, indices = indices)

}



pairs.transport <- function(Xfrom, Xto, eps, prefix){
  library(mop)


  gmraFrom <- c()
  indFrom = c()
  for(i in 1:length(Xfrom)){
    if(nrow(Xfrom[[i]]) > 0){
      
      gmraFrom <- c(gmraFrom, multiscale.transport.create.ipca(X=Xfrom[[i]], d=3,
          eps=eps, t=0.9) )
      indFrom = c(indFrom, i)
    }
  }

  gmraTo <- c()
  indTo = c()
  for(i in 1:length(Xto)){
    if(nrow(Xto[[i]]) > 0){
      
      gmraTo <- c(gmraTo, multiscale.transport.create.ipca(X=Xto[[i]], d=3,
          eps=eps, t=0.9) )
      indTo = c(indTo, i)
    }
  }

  
  dist = matrix( 0, nrow=length(gmraFrom), ncol=length(gmraTo) )

  for(i in 1:length(gmraFrom)){
    for(j in 1:length(gmraTo)){
      trp = multiscale.transport.id(gmra1=gmraFrom[i], gmra2 = gmraTo[j],
          p=2, rFactor=1, sType=0, scale1=-1, scale2=-1, stpPct=-1,
          oType=26, propFactor=0,  nRefinementIterations = 1 )
     
      dist[i, j] = trp$cost[length(trp$cost)]
#save(trp, file=sprintf("%s-%d-%d.Rdata", prefix, i, j))
    }
  }

  list(D = dist, indFrom = indFrom, indTo = indTo)

}




pairwise.assignments.segs <- function(segs){
 
  library(pdist) 
  library(lpSolve)

  D = matrix(0, nrow(segs), nrow(segs))
  for(i in 1:(nrow(segs)-1) ){
    print(i)
    for(j in (i+1):(nrow(segs)-1) ){
       pd = pdist(seg[i, ], seg[j, ])
       a = lp.assign(pd)
       D[i,j] = a$objval
       D[j,i] = a$objval
    }
  }

  D
}
