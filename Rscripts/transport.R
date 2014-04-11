##
#frame odor on:   5400
#frame odor off: 10801
#fly dimensions 25 x 12 pixels
#instantaneous velocity < 3 pixels /s , values abouve 4 p/s probably errorenous
#peak fluctation is 0.15 pixels, sd 0.09 pixels



#Pairwise Wassertsine distance for list of matrices in Xin
pairwise.transport <- function(Xin, eps=-1, scale=-1){
  library(mop)
  library(pdist)


  if( nrow(Xin[[1]]) > 200){

  gmra <- c()
  indices = c()
  for(i in 1:length(Xin)){
    nr = nrow(Xin[[i]])
    if(!is.null(nr)){
      if( nr > 0){

        gmra <- c(gmra, multiscale.transport.create.ipca(X=Xin[[i]], d=1,
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
    gc()
        print(i)
  }

  
    for(i in indices){
      multiscale.transport.delete.ipca(i)
    }

  res = list(D = dist)

  }
  else{
    dist = matrix( 0, nrow=length(Xin), ncol=length(Xin) )

      for(i in 1:(length(Xin)-1) ){
        for(j in (i+1):length(Xin)){
          C = pdist(Xin[[i]], Xin[[j]]);
          trp = transport(rep(1/nrow(Xin[[i]]), nrow(Xin[[i]])),
                rep(1/nrow(Xin[[j]]), nrow(Xin[[j]])), as.matrix(C) )
            dist[i, j] = trp$cost
            dist[j, i] = dist[i, j]
#save(trp, file=sprintf("%s-%d-%d.Rdata", prefix, i, j))
        }
        gc()
        print(i)
      }

  res = list(D = dist)

  }

  res
}



#Wasserstein distances between list of matrices Xfrom to Xto
pairs.transport <- function(Xfrom, Xto, eps=-1, scale=-1){
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
          p=2, rFactor=1, sType=0, scale1=scale, scale2=scale, stpPct=-1,
          oType=26, propFactor=0,  nRefinementIterations = 1 )
     
      dist[i, j] = trp$cost[length(trp$cost)]
#save(trp, file=sprintf("%s-%d-%d.Rdata", prefix, i, j))
    }
    gc()
  }

  list(D = dist, indFrom = indFrom, indTo = indTo)

}




#Check if the distance between X1 and X2 is statistically significant
transport.permutation.test <- function( X1, X2, eps, scale=0, nPerms=100 ){
  library(mop)

  trp = multiscale.transport.ipca(X1=X1, X2=X2, eps1=eps, eps2=eps, d1=3, d2=3,
      p=1, rFactor=1, sType=0, scale1=scale, scale2=scale, stpPct=-1, oType=26,
      propFactor=0,  nRefinementIterations = 1 )

  c1 = trp$cost
  
  X = rbind(X1, X2);
  n1 = nrow(X1)
  n2 = nrow(X2)
 
  nLess = 0

  cc = c();
  for(i in 1:nPerms){
    index= order(runif(n1+n2))
    Xa = X[index[1:n1], ]
    Xb = X[index[(n1+1):(n1+n2)], ]  
    
    trp = multiscale.transport.ipca(X1=Xa, X2=Xb, eps1=eps, eps2=eps, d1=3, d2=3,
      p=1, rFactor=1, sType=0, scale1=scale, scale2=scale, stpPct=-1, oType=26,
      propFactor=0,  nRefinementIterations = 1 )

     cc  = c(cc, trp$cost)

  }

  sum(cc > c1)/nPerms 

}



#
transport.pairwise.seg.len <- function( X, seg.len = c(10,20,40,80,160,320,640, 1280),
    seg.off = c(10,10,20,40,80,160,320, 640)  ){

  dists <- list();
  for(i in 1:length(seg.len) ){
     S = extract.markov.sample(X, seg.len[i], seg.off[i] )
     trp <- pairwise.transport(S$segs, -1, -1)
     dists[[i]] = trp$D 
  }

  dists
}


transport.pairwise.markov.len <- function( fa, index, ml = c(1,2,4,8), seg.len =
    50, seg.off = 25  ){

  dists <- list();
  for(i in 1:length(ml) ){
     tp = extract.markov.transitions(fa, ml[i])
     S = extract.markov.sample(tp[[index]], seg.len, seg.off )
     trp <- pairwise.transport(S$segs, -1, -1)
     dists[[i]] = trp$D 
  }

  dists
}




transport.agglomerate <- function( X, pCut = 0.05, eps=0, scale=0, nPerms=250){
  
  startIndex = 1
  compareIndex = 2

  n = length(X)

  agg = c()
  aggIndex = 1;
  
  X1 = X[[startIndex]]
  while(compareIndex <= n){
    p = transport.permutation.test(X1, X[[compareIndex]], eps, scale, nPerms);

    if( p > pCut ){
      X1 = rbind(X1, X[[compareIndex]])
      compareIndex = compareIndex + 1

    }
    else if(compareIndex < n){
      agg = rbind( agg, c(startIndex, compareIndex-1) )
      startIndex = compareIndex
      compareIndex = compareIndex + 1
      X1 = X[[startIndex]]
    }
    else{
      agg = rbind( agg, c(n, n) )
      compareIndex = compareIndex+1
    } 

  
  }


  agg

}

transport.agglomerate.iterative <- function ( S, nIter = 10, pCut = 0.05, eps=0, scale=0, nPerms=250){
    agg <- list()
    Sagg <- list()

    Sprev = S;
    for(i in 1:nIter){
      agg[[i]] = transport.agglomerate(Sprev$segs, pCut, eps, scale, nPerms)
      Sagg[[i]] = combine.agglomerate.segments(Sprev, agg[[i]])
      Sprev= Sagg[[i]]
    }
    A = list(Sagg=Sagg, agg=agg)
    A
}



combine.agglomerate.segments <- function(S, agg){
  
  conds = list()
  segs = list()
  for(i in 1:nrow(agg) ){
    segs[[i]] = S$segs[[agg[i,1]]]
    conds[[i]] = S$conds[[agg[i,1]]]
    if((agg[i,2] - agg[i, 1]) > 0){
      for(j in (agg[i,1]+1):agg[i,2]){
        segs[[i]] = rbind(segs[[i]], S$segs[[j]])
        conds[[i]] = rbind(conds[[i]], S$conds[[j]])
      }
    }
  }

  Sagg = list(segs=segs, conds=conds)

    Sagg
}






#cluster list of matrices in X based on statistical signficance of their
#Wasserstein distance
transport.cluster <- function(X, pCut = 0.01, nPerms=200, eps=0, scale=0){
  clusters <- list()
  cIds <- list();
  cIds = rep(-1, length(X));

  cc=1

  while(min(cIds) < 0){
    index = which(cIds < 0)[1]
    clusters[[cc]] = X[[index]]
    cIds[index] = cc

    for(i in 1:length(X)){
      if(cIds[i] < 0){
        p = transport.permutation.test(clusters[[cc]], X[[i]], eps, scale, nPerms)
        if(p > pCut){
          cIds[i] = cc
          clusters[[cc]] = rbind(clusters[[cc]], X[[i]])
        }
      }
    }

    cc = cc + 1;
  }

  list(clusters=clusters, cIds = cIds)
}
