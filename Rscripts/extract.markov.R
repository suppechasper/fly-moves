
#expects a list of data frames extract with extract.features
#k length of markov memory
extract.markov.transitions <- function(data, k){
  M = list()
  for( i in 1:length(data) ){
    n = length( data[[i]]$lengths )
    for(j in 1:k){
      if(j==1){

#tp= data[[i]]$ris[j:(n-k+j), ] 
tp = cbind(data[[i]]$lengths[j:(n-k+j)],
            data[[i]]$curvatureMean[j:(n-k+j)] )
      }
      else{
#        tp = cbind(tp,  data[[i]]$ris[ ,j:(n-k+j)] )
tp = cbind(tp,  data[[i]]$lengths[j:(n-k+j)], data[[i]]$curvatureMean[j:(n-k+j)] )
      }
    }

    dirim = 0
    dorim = 0
    time= 0
    x = 0
    for(j in 1:k){
     dirim = dirim + data[[i]]$dirim[j:(n-k+j)]
     dorim = dorim + data[[i]]$dorim[j:(n-k+j)]
     time = time + data[[i]]$time[j:(n-k+j)]
     x = x + data[[i]]$x[j:(n-k+j), ]
    }
    dirim = dirim/k
    dorim = dorim/k
    time = time/k
    x = x/k


    conds = cbind(  dirim, dorim, time, x, i) 
    colnames(conds) <-  c("dirim", "dorim", "time", "x1", "x2", "id") 
    

    #remove transition probabilites that contain a jump due to missing data
    span = data[[i]]$time[k:n] - data[[i]]$time[1:(n-k+1)]

    index = which(span < k)
    conds = conds[index, ]
    tp = tp[index, ]

    M[[i]] = list(tp = tp, conds= conds)
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





