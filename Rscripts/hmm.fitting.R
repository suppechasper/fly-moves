fit.hmm <- function(K, S, nStates, rtpm=F, rRho=F){

  library(hmm.discnp)

  Ktotal = c()

  for(i in 1:length(S$conds)){
    Ktotal = c(Ktotal, rep(K[i], nrow(S$conds[[i]]) ) )
  }

  mm = hmm(Ktotal, K=nStates, itmax=1000, rand.start=list(tpm=rtpm, Rho=rRho))
  print( mm$log.like )

  seq = mps(object = mm)


  seq


}



fit.depmix <- function(X, nstates){
  library(depmix)

  dep <- dmm(nstates=nstates, itemtypes = rep(1, ncol(X)) )
  md = markovdata(X, itemtypes = rep("continuous", ncol(X)) )

  depm = fitdmm(md, dep, printlevel=100);

}
