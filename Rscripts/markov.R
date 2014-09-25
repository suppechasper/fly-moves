normalize.rows <- function(M){
 
  for(i in 1:nrow(M)){
    M[i, ] = M[i, ] / sum(M[i, ])
  }
  M[is.na(M)]  = 0

  M

}

#build markov transition matrix from a sequence of cluster ids (state) and times
build.markov.transition <- function(clusters, time, nstates=max(clusters) ){

  M = matrix(0, nrow=nstates, ncol=nstates)
  for(i in 2:length(clusters)){
    if(time[i]-time[i-1] < 2){
      M[clusters[i-1], clusters[i] ]  = M[clusters[i-1], clusters[i] ] + 1
    }
  }
 
  normalize.rows(M)
}


#generate a path of length length from the transition matrix P with start state state 
#each state is mapped to a specific speed and curvature from which the path is
#generated accoridng to the transitions in P

generate.path <- function(P, state, length, speed, curvature){
  X <- matrix(0, nrow=length+1, ncol=2)
  X[2, ] = X[1, ] + speed[state]
  for(i in 2:length){
    cs = cumsum( P[state, ] )
    s = runif(1)
    for(j in 1:length(cs)){
      if(s < cs[j]){
        state = j;
        break;
      }
    }
    v = X[i, ] - X[i-1, ]
    phi = curvature[state]
    R = matrix( c(cos(phi), sin(phi), -sin(phi), cos(phi)), nrow=2)
    v = v / sqrt(sum(v^2))
    v = R %*% v
    X[i+1, ] = X[i, ]+v*speed[state] 
  }

  X
}



