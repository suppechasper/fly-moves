
extract.all.density.samples <- function(tp, len, off=len/3, n=10){

  D = c()
  conds = c()

  for(i in 1:length(tp)){
    S = extract.density.sample(tp[[i]], len, off, n);

    D = rbind(D, S$D)

    C = matrix(tp[[i]]$conds$id[1], nrow=nrow(S$D), ncol=12)
    
    start1 = which(tp[[i]]$conds$time > 5400 & tp[[i]]$conds$dirim > 0)
    start1 = tp[[i]]$conds$time[start1[1]]
    start2 = start1+300
    start3 = start2+600
    for( k in 1:length(S$conds) ){
      C[k, 2] = sum( S$conds[[k]]$time < start1 )
      C[k, 3] = sum( S$conds[[k]]$time > start1 & S$conds[[k]]$time < start2 )
      C[k, 4] = sum( S$conds[[k]]$time > start2 & S$conds[[k]]$time < start3 )
      C[k, 5] = sum( S$conds[[k]]$time > start3 & S$conds[[k]]$time < 10800 )
      C[k, 6] = sum( S$conds[[k]]$time > 10800 )
      C[k, 7] = which.max(C[k, 2:6])
      C[k, 8] = mean(S$conds[[k]]$time) - start1
      
     C[k, 9] =  sum(S$conds[[k]]$dorim < 20)
     C[k, 10] = sum( S$conds[[k]]$dorim >= 20 & S$conds[[k]]$dirim < -5)
     C[k, 11] = sum( S$conds[[k]]$dirim > -5 )
      C[k, 12] = which.max(C[k, 9:11])
    }
    conds = rbind(conds, C)
    
  }

  list(D=D, conds=conds)
}




extract.density.sample <- function( M, len, off=len/3, n=10 ){
  library(MASS)

  d = dim(M$tp)
  segs = list()
  conds = list()
  index = 1
  s = seq(1,(d[1]-len-1), off)
  D = matrix(0, nrow=length(s), ncol=n^2)
  for( i in s ){
    X = M$tp[i:(i+len-1), ]
    p = kde2d(X[,1], X[,2], h=1, lims=c(0, 6,-3.2, 3.2),  n=n  ) 

    segs[[index]] = p
    conds[[index]] = M$conds[i:(i+len-1), ]
    D[index, ] = as.vector(p[[3]])
    
    index = index +1  
  }
  
  res = list(segs=segs, conds=conds, D=D)

 res

}



