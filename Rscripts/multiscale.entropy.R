
extract.segs <- function(X, time, k, offset = k/2, rm.jumps=TRUE){
  

  n = nrow(X)
  S1 = c()
  S2 = c()
  for(i in 1:ncol(X)){
    for(j in 1:(k+offset)){
      s = seq(j, n-k-offset+j, by=offset)
      S2 = cbind( S2, X[s, i])
      if(j <= k){
         S1 = cbind( S1, X[s, i])
      }
    }
  } 
  if(rm.jumps){
    s = seq(1, n-k-offset+1, by=offset)
    tj = time[s+k+offset] - time[s] 
    ind  = which(tj == ( k+offset))
    S1 = S1[ind, ]
    S2 = S2[ind, ]
  }

  list(S1=S1, S2=S2)
}


multiscale.entropy <- function(X, time, nScales=5, r=0.15, threshold = -1){

  source("extract.features.R")

  sigma=1
  Xs = X
  ent = rep(NA, nScales)
  for(i in 1:nScales){
    S = extract.segs(Xs, time, k= 2^i )
    D1 <- dist(S$S1, method="maximum")
    D2 <- dist(S$S2, method="maximum")
    if(threshold < 0){
      threshold = sqrt(var(D1))*r
    }
    C1 = rowSums( as.matrix(D1) < threshold  ) - 1
    C2 = rowSums( as.matrix(D2) < threshold  ) - 1
    ent[i] = log( sum(C1) / sum(C2) )
    for(j in 1:ncol(X) ){
      Xs[,j] = ksmooth(time, X[,j], kernel="normal", bandwidth=2^(i-1)*4 )$y
    }
  }
  
  ent
}

