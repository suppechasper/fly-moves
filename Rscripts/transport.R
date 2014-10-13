#--- Pairwise optimal transport ---#


#Pairwise Wasserstein distance for list of matrices in Xin
pairwise.transport <- function(Xin, eps=-1, scale=-1, d=2, store.plan = FALSE,
    p=2, lambda=0, oType=26, rFactor=2, stop=6, split=1, sType=0,
    weight=as.list(rep(1, length(Xin))) ){


  library( gmra )
  library( mop )
  library( pdist )


  if( nrow(Xin[[1]]) > 200){

  gmra <- list()
 
  indices = c()
  for(i in 1:length(Xin)){
    nr = nrow(Xin[[i]])
    if(!is.null(nr)){
      if( nr > 0){

        gmra <- c(gmra, gmra.create.ipca(X=Xin[[i]], d=d,
            eps=eps, t=0.9, split=split, stop=stop) )
        indices = c(indices, i)
      }
    }
  }
 
  
  
 
  plans = list() 
  dist = matrix( 0, nrow=length(gmra), ncol=length(gmra) )

  count = 1;
  for(i in 1:(length(gmra)-1)){
    for(j in (i+1):length(gmra)){

      trp = multiscale.transport(gmra1=gmra[[i]], gmra2 = gmra[[j]], p=p,
          scale1=scale, scale2=scale, matchScale=FALSE, rFactor=rFactor,
          sType=sType, lambda=lambda, oType=oType, w1=weight[[i]], w2=weight[[j]])
        
#      trp = multiscale.transport.randomized.id(gmra1=gmra[i], gmra2 = gmra[j], p=p,
#          scale1=scale, scale2=scale, oType=26, matchScale=FALSE, nTrials=5)
     
      dist[i, j] = trp$cost[length(trp$cost)]
      dist[j, i] = dist[i, j]

      if(store.plan){
        plans[[count]] = list(trp = trp, i=i, j=j)
      }
      count=count+1
      print(paste(i, j))
#save(trp, file=sprintf("%s-%d-%d.Rdata", prefix, i, j))
    }
    gc()
  }

  
  res = list(D = dist, plans= plans)

  }
  else{
    plans = list() 
    dist = matrix( 0, nrow=length(Xin), ncol=length(Xin) )

      for(i in 1:(length(Xin)-1) ){
        for(j in (i+1):length(Xin)){
          C = as.matrix(dist(Xin[[i]], Xin[[j]]))^p
          trp = transport(rep(1/nrow(Xin[[i]]), nrow(Xin[[i]])),
                rep(1/nrow(Xin[[j]]), nrow(Xin[[j]])), as.matrix(C), lambda=lambda, oType=oType )
            dist[i, j] = trp$cost
            dist[j, i] = dist[i, j]    
            if(store.plan){
              plans[[count]] = list(trp = trp, i=i, j=j)
            }
            count=count+1

#save(trp, file=sprintf("%s-%d-%d.Rdata", prefix, i, j))
        }
        magc()
        print(i)
      }

  res = list(D = dist, plans =plans)

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



transport.extract.paths <- function(trp, scale = length(trp$cost) ){
  plan = trp$map[[scale]]

  paths = trp$to[[scale]][plan[,2], ] -  trp$from[[scale]][plan[,1], ] 

  res <- list(paths=paths, weights = plan[,3])
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






pairwise.transport.time <- function(tp, scale, eps){

  
  X = list();

  index = 1;
  for(i in 1:length(tp)){
    index1 = tp[[i]]$conds$time < 5400
    index2 = tp[[i]]$conds$time > 5400 & tp[[i]]$conds$time < 10800
    index3 = tp[[i]]$conds$time > 10800

    X[[index]] = tp[[i]]$tp[index1, ]
    X[[index+1]] = tp[[i]]$tp[index2, ]
    X[[index+2]] = tp[[i]]$tp[index3, ]

    index = index+3;
  }


  trp = pairwise.transport(X, scale=scale, eps=eps)
  res = list(X = X, trp =trp)

  res
}


plot.pc.transport.map <-function( ot, index, pcs = c(1,2) ){

  plan = ot$map[[index]]
  V = ot$from[[index]][ plan[,1], ] - ot$to[[index]][ plan[,2], ]
  
  V = V# * rep(plan[,3], ncol(V))

  pc <- prcomp(V, center=FALSE)
  X1 = ot$from[[index]] %*% pc$rotation
  X2 = ot$to[[index]] %*% pc$rotation
    
  multiscale.transport.plot.map(ot, index, pointAlpha=0, arrows=T,
        mapAlpha=1, arrow.angle=10, arrow.length=0.1, lwd=2,
        xlab= sprintf("OT-PC %d", pcs[1]), ylab= sprintf("OT-PC %d", pcs[2]), X1=X1,
        X2=X2 )

  pc
}




plot.binnned.transport.map <- function(trp, xbins, ybins, xlab="Curvature",
    ylab="Speed", col="black", arrow.length=0.1, arrow.angle=15, alpha=1, lwd=2,
    add = FALSE, cex.axis=1, cex.lab=1, cex=1, useTransparancy=FALSE,
    useCost=FALSE, maxW, maxC){

 
  l = length(trp$cost)


  nx = length(xbins)
  xloc = xbins[1:(nx-1)] + (xbins[2:nx] - xbins[1:(nx-1)]) / 2 
  
  ny = length(ybins)
  yloc = ybins[1:(ny-1)] + (ybins[2:ny] - ybins[1:(ny-1)]) / 2 


  w = matrix(0, ncol = length(xbins)-1, nrow = length(ybins)-1)
  cost = w;
  dx = w
  dy = w

  plan = trp$map[[l]]
  from = trp$from[[l]]
  to = trp$to[[l]]
  
  xb = cut(from[ plan[, 1], 1] , xbins)
  yb = cut(from[ plan[, 1], 2], ybins)

  delta = to[ plan[, 2], ] - from[ plan[,1], ]
  for(i in 1:nrow(plan)){
    cost[yb[i], xb[i]] = cost[yb[i], xb[i]]  + plan[i, 3] * plan[i, 4]
      
    w[yb[i], xb[i]] = w[yb[i], xb[i]]  + plan[i, 3]
    dx[yb[i], xb[i]] = dx[yb[i], xb[i]]  + plan[i, 3] * delta[i, 1]
    dy[yb[i], xb[i]] = dy[yb[i], xb[i]]  + plan[i, 3] * delta[i, 2]
  }

  dx = dx / w
  dy = dy / w
  w = w/max(w)
  cost = cost/max(cost)

  if(useCost){
    w=cost
  }

  if(!add){
    plot(NA, xlim=range(xbins), ylim=range(ybins), xlab=xlab, ylab=ylab,
        cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, bty="n")
  }

  for(i in 1:ncol(w)){
    for(j in 1:nrow(w)){
      if(useTransparancy){
        arrows( x0 = xloc[i], y0 = yloc[j], x1 = xloc[i]+dx[j, i], y1 = yloc[j]
          +dy[j, i], col=rgb( t(col2rgb(col)/255), alpha=0.05+0.95*w[j,i]), lwd=lwd, angle=arrow.angle,
          code=2, length=arrow.length )
      }
      else{
        arrows( x0 = xloc[i], y0 = yloc[j], x1 = xloc[i]+dx[j, i], y1 = yloc[j]
          +dy[j, i], col=rgb( t(col2rgb(col)/255), alpha=0.7), lwd=1+w[j,i]*lwd, angle=arrow.angle,
          code=2, length=arrow.length )
      }
    }
  }

  list(w=w, dx=dx, dy=dy)

}
