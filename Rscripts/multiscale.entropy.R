
extract.segs <- function(X, time, k){
  
  n = nrows(X)
  S = c()
  for(i in 1:ncol(X)){
    for(j in 1:k){
      s = seq(j, n-k+j, by=offset)
      S = cbind( S, X$[s, i])
    }
  } 
  if(rm.jumps){
    tj = time[k:n] - time[1:(n-k+1)] 
    s = seq(j, n-k+j, by=offset)
  
    ind  = which(jump == ( k-1))
    S = S[ind, ]
  }

  S
}

multiscale.entropy <- function(X, time, max.Smooth=16){

  source("extract.features.R")

  threshold


  sigma=1
  Xs = X
  for(i in 1:max.Smooth){
    extract.segs(X
    D <- dist(X, method="Euclidean")


}

