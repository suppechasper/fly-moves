
hyperbolic.mds <- function(D, m, k){

  C = cosh(sqrt(k)*D)
  eig = eigen(C, TRUE)

  d = length(eig$values)
  index = c((d-m+1):d)

  A = eig$vector[,index] %*% diag( sqrt(abs(eig$values[index])))

# s = sqrt( A[,1]^2 - rowSums(A[,2:(m+1)]^2) )
# x = sweep(A,MARGIN=1,s,'/')
# browser()
  hmds = list( x = cbind(1+rowSums(A^2), A), eig=eig)
}






spherical.mds <- function(D, m, k){

  C = cos(sqrt(k)*D)
  eig = eigen(C, TRUE)

  d = length(eig$values)
  index = c(1:(m+1))

  A = eig$vector[,index] %*% diag(sqrt(abs(eig$values[index])))

# s = rowSums(A^2)
# x = sweep(A, MARGIN=1, s,'/')
  smds = list( x=A, eig=eig)
}







spherical.distance <- function(X, k){
 D=matrix(0, nrow=nrow(X), ncol=nrow(X))
 for(i in 1:nrow(X)){
   D[,i] = X %*% X[i,]
 }

 D = acos(D)/sqrt(k)

}






minkowski.mds <- function(D, p, n=ncol(X)-p){

  D2 = D^2
  #double centering
  C = sweep( D2, 1, rowMeans(D2), "-" )
  C = -0.5 * sweep( C, 2, colMeans(C), "-" )

  eig = eigen(C, TRUE)

  d = nrow(C)
  index = c(1:p, (d-n+1):d)

  A = eig$vector[,index] %*% diag( sqrt(abs(eig$values[index]) ) )

  list(x=A, eig=eig$values)

}






#hyperboloid distance for metric (1, -1, ..., -1)
minkowski.distance <- function(X, p, n=ncol(X)-p){ 
  
 M = diag(c(rep(1, p), rep(-1, n) ) )

 D=matrix(0, nrow=nrow(X), ncol=nrow(X))
 for(i in 1:nrow(X)){
#Xc = t(apply(X, 1, '-', X[i,]))
#  D[,i] = rowSums( (Xc %*% M ) * Xc ) 
   D[,i] = X %*% M  %*% X[i, ] 
 }

 D 

}





#hyperboloid distance for metric (1, -1, ..., -1)
hyperboloid.distance <- function(X, k){ 
  
 M = diag(rep(-1, ncol(X)))
 M[1,1]=1


 D=matrix(0, nrow=nrow(X), ncol=nrow(X))
 for(i in 1:nrow(X)){
   D[,i] = X %*% M %*% X[i,]
 }

 browser()
 D =  acosh(D)/sqrt(k)

}





#pairwse distance from hyperbolic data
generate.hyperboloid.data <-function(m, n, s=0.1){
 X = matrix(runif(m*n)-0.5, nrow=m, ncol=n)*s
 
 X= cbind( 1+rowSums(X^2), X )
 
 D =  hyperboloid.distance(X)
}
