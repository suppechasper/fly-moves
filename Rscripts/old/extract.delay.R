
#expects a list of data frames extract with extract.features
#k length of markov memory
extract.markov.transitions <- function(X, k){
    n = length( X$lengths )
    for(j in 1:k){
      if(j==1){

#tp= data[[i]]$ris[j:(n-k+j), ] 
tp = cbind(X$lengths[j:(n-k+j)],
            X$curvatureMean[j:(n-k+j)] )
      }
      else{
#        tp = cbind(tp,  data[[i]]$ris[ ,j:(n-k+j)] )
tp = cbind(tp,  X$lengths[j:(n-k+j)], X$curvatureMean[j:(n-k+j)] )
      }
    }

    dirim = 0
    dorim = 0
    time= 0
    x = 0
    for(j in 1:k){
     dirim = dirim + X$dirim[j:(n-k+j)]
     dorim = dorim + X$dorim[j:(n-k+j)]
     time = time + X$time[j:(n-k+j)]
     x = x + X$x[j:(n-k+j), ]
    }
    dirim = dirim/k
    dorim = dorim/k
    time = time/k
    x = x/k


    conds = cbind(  dirim, dorim, time, x, i) 
    colnames(conds) <-  c("dirim", "dorim", "time", "x1", "x2", "id") 
    

    #remove transition probabilites that contain a jump due to missing data
    span = X$time[k:n] - X$time[1:(n-k+1)]

    index = which(span < k)
    conds = conds[index, ]
    tp = tp[index, ]

    list(tp = tp, conds= conds)
}




extract.all.delay <- function(data, k){
  M = list()
  for( i in 1:length(data) ){
    d = extract.delay(data[[i]], k)
    M[[i]] = d
  }

  M
}


#extract a sample of the markov transition probabilites
# M is a single elemnt of the list extract with extract.markov transition
#len is the size of the sample to extract
extract.markov.sample <- function( M, len, off=len/3 ){
  d = dim(M$tp)
  segs = list()
  conds = list()
  index = 1
  for( i in seq(1,(d[1]-len-1), off) ){
    segs[[index]] = M$tp[i:(i+len-1), ]
    conds[[index]] = M$conds[i:(i+len-1), ]
    index = index +1  
  }
  
  res = list(segs=segs, conds=conds)

 res

}

extract.all.markov.sample <- function( tp, len, off=len/3 ){
  segs = list()
  conds = list()
  index = 1
  for(j in 1:length(tp)){
    M = tp[[j]]
  
    d = dim(M$tp)
    for( i in seq(1,(d[1]-len-1), off) ){
      segs[[index]] = M$tp[i:(i+len-1), ]
      conds[[index]] = M$conds[i:(i+len-1), ]
      index = index +1  
    }
  
  }
  res = list(segs=segs, conds=conds)

 res

}

