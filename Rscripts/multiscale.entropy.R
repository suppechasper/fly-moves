
extract.segs <- function(X, time, k){
  
  n = nrows(X)
  S1 = c()
  S2 = c()
  for(i in 1:ncol(X)){
    for(j in 1:k){
      s = seq(j, n-k+j, by=k)
      S1 = cbind( S1, X$[s, i])
      if(j != k){
        S2 = cbind( S2, X$[s, i])
      }
    }
  } 
  if(rm.jumps){
    tj = time[k:n] - time[1:(n-k+1)] 
    ind = which(jump == ( k-1))
    S1 = S1[ind, ]
    S2 = S2[ind, ]
  }

  list(S1=S1, S2=S2)
}

multiscale.entropy <- function(X, k=3,time, max.Smooth=16){

  source("extract.features.R")

  threshold = 0;
  sigma=1
  Xs = X
  for(i in 1:max.Smooth){
    S = extract.segs(X, time, k)
    D1 <- dist(S$S1, method="Maximum")
    D2 <- dist(S$S2, method="Maximum")
    if(threshold = 0){
      threshold = sqrt( var(D1) ) * 0.15
    }
    
    C1 = rowSums(as.matrix(D1) < threshold

}

