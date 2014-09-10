
delay.reconstruction <- function(X, k){
  
    n = nrow( X )
    delay = c();
    for(i in 1:ncol(X)){
      for(j in 1:k){
        delay = cbind(delay,  X[j:(n-k+j), i]) 
      }
    }

    delay
}
