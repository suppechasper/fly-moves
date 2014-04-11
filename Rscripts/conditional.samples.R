##
#frame odor on:   5400
#frame odor off: 10801
#fly dimensions 25 x 12 pixels
#instantaneous velocity < 3 pixels /s , values abouve 4 p/s probably errorenous
#peak fluctation is 0.15 pixels, sd 0.09 pixels

conditional.samples <- function(data, conds, n=30){

  odor = conds$time > 5400 & conds$time < 10801
  inside = conds$dirim > -25
  atRim = abs(conds$dorim) < 25
  id = conds$id 


  #build  a set of conditional samples
  X <- list()
  
  C <- matrix(0, nrow= n * 18 + 18, ncol=6)
  colnames(C) <- c("odorOn", "odorOff", "innerRim", "middle", "outerRim", "id") 


  index = 1;
  #conditional averages over all flies
  X[[index]] = data
  C[index, ] = 1
  index = index + 1
  
  #
  X[[index]] = data[inside == 1, ]
  C[index, ] = c(1, 1, 1, 0, 0, -1)
  index = index +1
  
  X[[index]] = data[inside == 0 & atRim == 0, ]
  C[index, ] = c(1, 1, 0, 1, 0, -1)
  index = index +1

  X[[index]] = data[ inside == 0, ]
  C[index, ] = c(1, 1, 0, 1, 1, -1)
  index = index +1
  
  X[[index]] = data[ atRim==1, ]
  C[index, ] = c(1, 1, 0, 0, 1, -1)
  index = index + 1

  X[[index]] = data[ atRim==0, ]
  C[index, ] = c(1, 1, 1, 1, 0, -1)
  index = index + 1


  #
  X[[index]] = data[odor == 1 & inside == 1, ]
  C[index, ] = c(1, 0, 1, 0, 0, -1)
  index = index +1
  
  X[[index]] = data[odor == 1 & inside == 0 & atRim == 0, ]
  C[index, ] = c(1, 0, 0, 1, 0, -1)
  index = index +1

  X[[index]] = data[odor == 1 & inside == 0, ]
  C[index, ] = c(1, 0, 0, 1, 1, -1)
  index = index +1
  
  X[[index]] = data[odor==1 & atRim==1, ]
  C[index, ] = c(1, 0, 0, 0, 1, -1)
  index = index +1

  X[[index]] = data[odor==1, ]
  C[index, ] = c(1, 0, 1, 1, 1, -1)
  index = index +1

  X[[index]] = data[ odor==1 & atRim==0, ]
  C[index, ] = c(1, 0, 1, 1, 0, -1)
  index = index + 1


  #
  X[[index]] = data[odor == 0 & inside == 1, ]
  C[index, ] = c(0, 1, 1, 0, 0, -1)
  index = index +1
  
  X[[index]] = data[odor == 0 & inside == 0 & atRim == 0, ]
  C[index, ] = c(0, 1, 0, 1, 0, -1)
  index = index +1

  X[[index]] = data[odor == 0 & inside == 0, ]
  C[index, ] = c(0, 1, 0, 1, 1, -1)
  index = index +1
  
  X[[index]] = data[odor==0 & atRim==1, ]
  C[index, ] = c(0, 1, 0, 0, 1, -1)
  index = index +1

  X[[index]] = data[odor==0, ]
  C[index, ] = c(0, 1, 1, 1, 1, -1)
  index = index +1

  X[[index]] = data[ odor==0 & atRim==0, ]
  C[index, ] = c(0, 1, 1, 1, 0, -1)
  index = index + 1






  #per fly conditions
  for(i in 1:max(conds$id) ){
    print(i)

  X[[index]] = data[id == i, ]
  C[index, ] = c(1, 1, 1, 1, 1, i)
  index = index + 1
  
  #
  X[[index]] = data[inside == 1 & id == i, ]
  C[index, ] = c(1, 1, 1, 0, 0, i)
  index = index +1
  
  X[[index]] = data[inside == 0 & atRim == 0 & id == i, ]
  C[index, ] = c(1, 1, 0, 1, 0, i)
  index = index +1

  X[[index]] = data[ inside == 0 & id == i, ]
  C[index, ] = c(1, 1, 0, 1, 1, i)
  index = index +1
  
  X[[index]] = data[atRim==1 & id == i, ]
  C[index, ] = c(1, 1, 0, 0, 1, i)
  index = index +1
  
  X[[index]] = data[ atRim==0 & id==i, ]
  C[index, ] = c(1, 1, 1, 1, 0, i)
  index = index + 1


  #
  X[[index]] = data[odor == 1 & inside == 1 & id == i, ]
  C[index, ] = c(1, 0, 1, 0, 0, i)
  index = index +1
  
  X[[index]] = data[odor == 1 & inside == 0 & atRim == 0 & id == i, ]
  C[index, ] = c(1, 0, 0, 1, 0, i)
  index = index +1

  X[[index]] = data[odor == 1 & inside == 0 & id == i, ]
  C[index, ] = c(1, 0, 0, 1, 1, i)
  index = index +1
  
  X[[index]] = data[odor==1 & atRim==1 & id == i, ]
  C[index, ] = c(1, 0, 0, 0, 1, i)
  index = index +1

  X[[index]] = data[odor==1 & id == i, ]
  C[index, ] = c(1, 0, 1, 1, 1, i)
  index = index +1

  X[[index]] = data[ odor==1 & atRim==0 & id==i, ]
  C[index, ] = c(1, 0, 1, 1, 0, i)
  index = index + 1


  #
  X[[index]] = data[odor == 0 & inside == 1 & id == i, ]
  C[index, ] = c(0, 1, 1, 0, 0, i)
  index = index +1
  
  X[[index]] = data[odor == 0 & inside == 0 & atRim == 0 & id == i, ]
  C[index, ] = c(0, 1, 0, 1, 0, i)
  index = index +1

  X[[index]] = data[odor == 0 & inside == 0 & id == i, ]
  C[index, ] = c(0, 1, 0, 1, 1, i)
  index = index +1
  
  X[[index]] = data[odor==0 & atRim==1 & id == i, ]
  C[index, ] = c(0, 1, 0, 0, 1, i)
  index = index +1

  X[[index]] = data[odor==0 & id == i, ]
  C[index, ] = c(0, 1, 1, 1, 1, i)
  index = index +1

  X[[index]] = data[ odor==0 & atRim==0 & id==i, ]
  C[index, ] = c(0, 1, 1, 1, 0, i)
  index = index + 1


  }

  n <- unlist( lapply(X, nrow) )

  res = list( samples = X, conditions=C, nSamples = n)

}




