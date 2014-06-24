#----- Extract data for analysis from raw data----


#---- Experimental setup information ----
#frame odor on:   5400
#frame odor off: 10801
#fly dimensions 25 x 12 pixels
#instantaneous velocity < 3 pixels /s , values abouve 4 p/s probably errorenous
#peak fluctation is 0.15 pixels, sd 0.09 pixels






#read file number index in above files list and extract features
extract.features <- function(xyFile, irFile, orFile,  index, lT=6, rm.na=T, std=0){

  x <- read.table( xyFile, sep=",", header=T)


  x = x + matrix(rnorm(nrow(x)*ncol(x), sd=std), ncol = ncol(x))

  irim = get.inner.rim(irFile);
  orim = get.outer.rim(orFile);
  n = nrow(x)

  v1 = x;
  v1[,1] = v1[,1] - irim$mean[1] 
  v1[,2] = v1[,2] - irim$mean[2] 
  dirim = irim$r - sqrt( rowSums(v1^2) )


  v2 = x;
  v2[,1] = v2[,1] - orim$mean[1] 
  v2[,2] = v2[,2] - orim$mean[2] 
  dorim = orim$r - sqrt( rowSums(v2^2) )

  v = x[2:n, ] - x[1:(n-1), ]
  vn = n-1
  step = sqrt( rowSums(v^2) )
  lengths = (step[1:(vn-1)] + step[2:vn] )/2


  angles = rowSums( v[1:(vn-1), ] * v[2:vn, ]) / (step[1:(vn-1)] * step[2:vn])
  angles[angles > 1] = 1
  angles[angles < -1] = -1
  cross = v[1:(vn-1), 1] * v[2:vn, 2] -  v[1:(vn-1), 2] * v[2:vn, 1]
  
  ris = cbind( angles, sign(cross)*sqrt(1-angles^2) ) * step[2:vn]; 

  time = 2:vn
  dorim = dorim[2:vn]
  dirim = dirim[2:vn]
  
  u = apply(cbind(step[1:(vn-1)], step[2:vn]), 1, max)
  if(rm.na == T){
  ind = which(u < lT)
  lengths = lengths[ind]
  angles  = angles[ind]
  cross   = cross[ind]
  time    = time[ind]
  dirim   = dirim[ind]
  dorim   = dorim[ind]
  x = x[ind +1, ]
  
  ind2 = complete.cases(angles)
  lengths = lengths[ind2]
  angles  = angles[ind2]
  cross   = cross[ind2]

  time    = time[ind2]
  dirim   = dirim[ind2]
  dorim   = dorim[ind2]
  x = x[ind2, ]
  }


  list( lengths = lengths, angles=angles, cross = cross, id = index, time = time,
        dorim = dorim, dirim = dirim, curvature = sign(cross)*acos(angles),
        ris=ris, x = x, u=u )

}


extract.all.features <- function(xyFiles, irFiles, orFiles, lT=6, rm.na=T){
  data <- list()
  for( i in 1:length(xyFiles) ){
    print(i)
    data[[i]] = extract.features(xyFiles[[i]], irFiles[[i]], orFiles[[i]], i, lT, rm.na)
  }
  data
}



#Extract features with noise added and avergae over nRuns
extract.all.features.expected <- function(xyFiles, irFiles, orFiles, lT=6, rm.na=T, nRuns=10, std=0.06){
  data <- list()
  for( i in 1:length(xyFiles)){
    print(i)
    dataSum = extract.features(xyFiles[[i]], irFiles[[i]], orFiles[[i]], i, lT, rm.na=F)
    dataSum$curvatureOrig = dataSum$curvature
    for(k in 1:nRuns){
      dataTmp = extract.features( xyFiles[[i]], irFiles[[i]], orFiles[[i]], i, lT, rm.na=F, std )
#data[[i]]$lengths = data[[i]]$lengths + dataTmp$lengths
      dataSum$cross = dataSum$cross + dataTmp$cross
      dataSum$u = dataSum$u + dataTmp$u
      dataSum$angles = dataSum$angles + dataTmp$angles
      dataSum$curvature = dataSum$curvature + dataTmp$curvature
      dataSum$ris = dataSum$ris + dataTmp$ris
    }
# data[[i]]$length = data[[i]]$length / (nRuns+1) 
    dataSum$cross = dataSum$cross / (nRuns+1) 
    dataSum$u = dataSum$u / (nRuns+1) 
    dataSum$angles = dataSum$angles / (nRuns+1) 
    dataSum$curvatureMean = dataSum$curvature / (nRuns+1) 
    dataSum$ris = dataSum$ris / (nRuns+1) 
  
    if(rm.na == T){
       ind = which(dataSum$u < lT)
        dataSum$lengths = dataSum$lengths[ind]
        dataSum$angles  = dataSum$angles[ind]
        dataSum$cross   = dataSum$cross[ind]
        dataSum$time    = dataSum$time[ind]
        dataSum$dirim   = dataSum$dirim[ind]
        dataSum$dorim   = dataSum$dorim[ind]
        dataSum$curvatureMean   = dataSum$curvatureMean[ind]
        dataSum$curvatureOrig   = dataSum$curvatureOrig[ind]
        dataSum$x = dataSum$x[ind +1, ]
        dataSum$ris  = dataSum$ris[ind, ]

        ind = complete.cases(dataSum$angles)
        dataSum$lengths = dataSum$lengths[ind]
        dataSum$angles  = dataSum$angles[ind]
        dataSum$cross   = dataSum$cross[ind]

        dataSum$time    = dataSum$time[ind]
        dataSum$dirim   = dataSum$dirim[ind]
        dataSum$dorim   = dataSum$dorim[ind]
        dataSum$curvatureMean   = dataSum$curvatureMean[ind]
        dataSum$curvatureOrig   = dataSum$curvatureOrig[ind]
        dataSum$x = dataSum$x[ind, ]
        dataSum$ris  = dataSum$ris[ind, ]

    }
    dataSum$curvature = sign(dataSum$cross)*acos(dataSum$angles)
    data[[i]] = dataSum

  }
  data
}




#get rim data for file i by circle fitting
get.inner.rim <- function(file){ 
  library(circular)
  x <- read.table( file, sep=",", header=F)
  
  x <- x[1:(nrow(x)-1), ]
  circ = lsfit.circle(x)
  rim = list( mean = circ$coefficients[2:3], r=circ$coefficients[1], x=x)

   rim

}



get.outer.rim <- function(file){
  library(circular)

  x <- read.table( file, sep=",", header=F)
  
  x <- x[1:(nrow(x)-1), ]
  circ = lsfit.circle(x)
  rim = list( mean = circ$coefficients[2:3], r=circ$coefficients[1], x=x)
  
  rim

}




#extract segments of length k from all features (read from the data with above
# methods)
extract.all.segments <- function(features, k){
  Slist = list()
  
  for(i in 1:length(features)){
    Slist[[i]] = extract.segments(features[[i]], k) 
    print(i)
  }

  Slist
}


#extract segements of length k from a single feature set X 
extract.segments  <- function(X, k){    
   n = length( X$lengths )
     lengths=c()
     curvature=c()
     dirim=c()
     dorim=c()
     time=c()
     x=c()
     y=c()
     rx=c()
     ry=c()
     for(j in 1:k){
        lengths = cbind( lengths, X$lengths[j:(n-k+j)])
        curvature = cbind( curvature, X$curvature[j:(n-k+j)])
        dirim = cbind( dirim, X$dirim[j:(n-k+j)])
        dorim = cbind( dorim, X$dorim[j:(n-k+j)])
        time = cbind( time, X$time[j:(n-k+j)])
        x= cbind( x,X$x[j:(n-k+j),1])
        y = cbind( y,X$x[j:(n-k+j),2])
        rx= cbind( rx,X$ris[j:(n-k+j),1])
        ry = cbind( ry,X$ris[j:(n-k+j),2])

    }
    C = matrix(0, nrow=nrow(lengths), ncol=11)

      start1 = which(X$time > 5400 & X$dirim > 0)
      start1 = X$time[start1[1]]
      start2 = start1+300
      start3 = start2+600
      for( k in 1:nrow(C) ){
        C[k, 1] = sum( X$conds[[k]]$time < start1 )
          C[k, 2] = sum( time[k, ] > start1 & time[k, ]  < start2 )
          C[k, 3] = sum( time[k, ] > start2 & time[k, ]  < start3 )
          C[k, 4] = sum( time[k, ] > start3 & time[k, ]  < 10800 )
          C[k, 5] = sum( time[k, ] > 10800 )
          C[k, 6] = which.max(C[k, 1:5])
          C[k, 7] = mean(time[k, ]) - start1

          C[k, 8] =  sum(dorim[k,] < 20)
          C[k, 9] = sum( dorim[k,] >= 20 & X$conds[[k]]$dirim < -5)
          C[k, 10] = sum( dirim[k,] > -5 )
          C[k, 11] = which.max(C[k, 8:10])
      }
    colnames(C) <- c("before", "odor1", "odor2", "odor3", "after", "maxOdor",
        "timeToOdor", "outer", "middler", "inner", "maxPos");

   list( lengths = lengths, curvature=curvature, dirim=dirim, dorim=dorim,
       time=time,x=x, y=y, rx=rx, ry=ry, C=as.data.frame(C), id=X$id)

}





extract.sample.segments <- function(S, len, off=len/2){
  d = dim(S$lengths)
  index = 1
  sample = list()
  for( i in seq(1,(d[1]-len-1), off) ){
    sample[[index]] = S[i:(i+len-1), ]
    index = index +1  
  }
  
 sample
} 





extract.sample <- function(S, len, off=len/2){
  d = dim(S$lengths)
  index = 1
  sample = list()
  for( i in seq(1,(d[1]-len-1), off) ){
    sample[[index]] = S[i:(i+len-1), ]
    index = index +1  
  }
  
 sample
} 



#given the data X with n rows and conditions C (i.e. the matrix C from
#extract.features) for each row extract submatrices condition on odor events
extract.condition.odor <- function(X, C){
  Xbefore = X[C$maxOdor == 1, ]
  nb = nrow(Xbefore)
  nbs = nb/2
  
  Xduring = X[C$maxOdor == 2 | C$maxOdor == 3 | C$maxOdor == 4, ]
  nd = nrow(Xduring)
  nds = nd/2

  Xafter = X[C$maxOdor == 5, ]
  na = nrow(Xafter)
  nas = na/2

  list( Xbefore1 = Xbefore[1:nbs, ], Xbefore2 = Xbefore[(nbs+1):nb, ],
        Xduring1 = Xduring[1:nds, ], Xduring2 = Xduring[(nds+1):nd, ],
        Xafter1 = Xafter[1:nas, ], Xafter2 = Xafter[(nas+1):na, ] )

}




extract.condition.odor.all <- function(Xlist, Slist){

  Xbefore1 = c()
  Xbefore2 = c()
  Xduring1 = c()
  Xduring2 = c()
  Xafter1 = c()
  Xafter2 = c()

  for(i in 1:length(Xlist) ){
    res = extract.condition.odor( Xlist[[i]], Slist[[i]]$C )
    Xbefore1 = rbind(Xbefore1, res$Xbefore1)
    Xbefore2 = rbind(Xbefore2, res$Xbefore2)
    Xduring1 = rbind(Xduring1, res$Xduring1)
    Xduring2 = rbind(Xduring2, res$Xduring2)
    Xafter1  = rbind(Xafter1, res$Xafter1)
    Xafter2  = rbind(Xafter2, res$Xafter2)
  }

  list(Xbefore1 = Xbefore1, Xbefore2 = Xbefore2, 
       Xduring1 = Xduring1, Xduring2 = Xduring2, 
       Xafter1 = Xafter1, Xafter2 = Xafter2)

}





extract.condition.odor.all.list <- function(Xlist, Slist){


  O <- list()
  index = 1
  for(i in 1:length(Xlist) ){
    res = extract.condition.odor( Xlist[[i]], Slist[[i]]$C )
    O[[index]] = res$Xbefore1;
    index = index +1
    O[[index]] = res$Xbefore2;
    index = index +1
    O[[index]] = res$Xduring1;
    index = index +1
    O[[index]] = res$Xduring2;
    index = index +1
    O[[index]] = res$Xafter1;
    index = index +1
    O[[index]] = res$Xafter2;
    index = index +1
  }

  O
}


