

complexity.stress.mds <- function(dists, dim=2){

  library(MASS)

  mds <- list()
  for(i in 1:length(dists) ){
    mds[[i]] <- isoMDS(dists[[i]], k=dim)
  }

  mds
}






complexity.clus.gap <- function(dists, K.max = 20, dim=3){
  
  library(MASS)
  library(cluster)

  cgap <- list()
  for(i in 1:length(dists) ){
    mds <- cmdscale(dists[[i]], k=dim)

    cgap[[i]] = clusGap(mds, kmeans, K.max = K.max, nstart=10, B=100)
  }

  cgap

}




complexity.clus.gap.get.clusters <- function(cgaps, factor=2, method="firstSEmax"){

  nc = c()
  for( i in 1:length(cgaps)){
    nc[i] = maxSE(cgaps[[i]]$Tab[,3], factor*cgaps[[i]]$Tab[,4], method=method)
  }

  nc
}






complexity.hclust.reconstruction.error <- function(dists, K.max=100,
    method="ward"){
  
  err <- matrix(0, nrow=K.max, ncol=length(dists) )

  for(i in 1:length(dists) ){
    hc <- hclust(as.dist( dists[[i]] ), method=method)
    for(k in 1:K.max){
      g = cut = cutree(hc, k=k)
      tmp = 0;
      for(j in 1:k){
        index = which(g==j)
        tmp =  tmp + sum(dists[[i]][index, index]) / length(index) 
      }
      err[k, i] = tmp
    }
  }

  err
}


auto.correlation <- function(dists, K.max = 100){

  A <- matrix(0, nrow=length(dists), ncol=K.max)
  for(n in 1:length(dists)) {
    m = mean(dists[[n]]) *  ncol(dists[[n]])^2 / ( ncol(dists[[n]])^2 + ncol(dists[[n]]) )

  for(k in 1:min( K.max, ncol(dists[[n]])-1 ) ){
    tmp = 0
    s = 1:(ncol(dists[[n]])-k)
    for(i in  s){
      tmp = tmp+dists[[n]][i, i+k]
    }
    A[n, k] = 1 - (tmp/length(s))/m
  }
  }

  A


}
