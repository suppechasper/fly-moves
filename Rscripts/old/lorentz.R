sample.lorentz <-function(xstart, n, a=10, b=28, c=8/3, t=0.01){

  x = xstart
  X = matrix(0, nrow=n, ncol=3);
  X[1,] = x
  for(i in 2:n){
    X[i,1] = a*(X[i-1,2] - X[i-1,1]) 
    X[i,2] = X[i-1,1]*(b - X[i-1,3]) - X[i-1,2] 
    X[i,3] = X[i-1,2]* X[i-1,1] - c*X[i-1,3]
    X[i,] = X[i-1,] + t*X[i,] 
  }

  X
}



#expects a list of data frames extract with extract.features
#k length of markov memory
delay.reconstruction <- function(x, k){
  
    n = length( x )
    delay = c();
    for(j in 1:k){
      delay = cbind(delay,  x[j:(n-k+j)] )
    }

    delay
}
