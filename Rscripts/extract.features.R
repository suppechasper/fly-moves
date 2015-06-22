#----- Extract data for analysis from raw data----


#---- Experimental setup information ----
#frame odor on:   5400
#frame odor off: 10801
#fly dimensions 25 x 12 pixels
#instantaneous velocity < 3 pixels /s , values abouve 4 p/s probably errorenous
#peak fluctation is 0.15 pixels, sd 0.09 pixels






#read file number index in above files list and extract features
extract.features <- function(xyFile, irFile, orFile,  index, lT=6, rm.na=T,
    std=0, uT=0){

  x <- read.table( xyFile, sep=",", header=T, colClasses="numeric")


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

  va = (v[1:(vn-1), ] + v[2:vn, ])/2.0
  dir = v2[1:(vn-1), ] + va
  orientation = acos( rowSums( va * dir) /( sqrt( rowSums(va^2)) * sqrt( rowSums(dir^2) ) ) )
  orientation[is.na(orientation)] = 0

  angles = rowSums( v[1:(vn-1), ] * v[2:vn, ]) / (step[1:(vn-1)] * step[2:vn])
  angles[angles > 1] = 1
  angles[angles < -1] = -1
  cross = v[1:(vn-1), 1] * v[2:vn, 2] -  v[1:(vn-1), 2] * v[2:vn, 1]
  
  ris = cbind( angles, sign(cross)*sqrt(1-angles^2) ) * step[2:vn]; 

  time = 2:vn
  dorim = dorim[2:vn]
  dirim = dirim[2:vn]
  
  u = apply(cbind(step[1:(vn-1)], step[2:vn]), 1, max)    
  v = apply(cbind(step[1:(vn-1)], step[2:vn]), 1, min)    
  

  x=x[2:(nrow(x)-1), ]

  if(rm.na == T){
    ind = which(u < lT)
      lengths = lengths[ind]
      angles  = angles[ind]
      cross   = cross[ind]
      time    = time[ind]
      dirim   = dirim[ind]
      dorim   = dorim[ind]
      ris = ris[ind, ]
      x = x[ind, ]
      orientation = orientation[ind]

      ind = which(v[ind] > uT)
      lengths = lengths[ind]
      angles  = angles[ind]
      cross   = cross[ind]
      time    = time[ind]
      dirim   = dirim[ind]
      dorim   = dorim[ind]
      ris = ris[ind, ]
      x = x[ind, ]
      orientation = orientation[ind]

      ind2 = complete.cases(angles)
      lengths = lengths[ind2]
      angles  = angles[ind2]
      cross   = cross[ind2]

      time    = time[ind2]
      dirim   = dirim[ind2]
      dorim   = dorim[ind2]
      x = x[ind2, ]
      ris = ris[ind2, ]
      orientation = orientation[ind2]
  }


  list( lengths = lengths, angles=angles, cross = cross, id = index, time = time,
        dorim = dorim, dirim = dirim, curvature = sign(cross)*acos(angles),
        ris=ris, x = x, u=u, v=v, irim=irim, orim=orim, orientation=orientation )

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
extract.all.features.expected <- function(xyFiles, irFiles, orFiles, lT=6,
    uT=0.5,  rm.na=T, nRuns=10, std=0.06){
  data <- list()
  for( i in 1:length(xyFiles)){
    print(i)
    dataSum = extract.features(xyFiles[[i]], irFiles[[i]], orFiles[[i]], i,
        lT=1000, rm.na=F, uT=-1)


    dataSum$curvatureOrig = dataSum$curvature
    for(k in 1:nRuns){
      dataTmp = extract.features( xyFiles[[i]], irFiles[[i]], orFiles[[i]], i,
          lT=1000, rm.na=F, std=std, uT=-1)
#data[[i]]$lengths = data[[i]]$lengths + dataTmp$lengths
        dataSum$cross = dataSum$cross + dataTmp$cross
        dataSum$u = dataSum$u + dataTmp$u
        dataSum$v = dataSum$v + dataTmp$v
        dataSum$angles = dataSum$angles + dataTmp$angles
        dataSum$curvature = dataSum$curvature + dataTmp$curvature
        dataSum$ris = dataSum$ris + dataTmp$ris
        dataSum$orientation = dataSum$orientation + dataTmp$orientation
    }
# data[[i]]$length = data[[i]]$length / (nRuns+1) 
    dataSum$cross = dataSum$cross / (nRuns+1) 
    dataSum$u = dataSum$u / (nRuns+1) 
    dataSum$v = dataSum$v / (nRuns+1) 
    dataSum$angles = dataSum$angles / (nRuns+1) 
    dataSum$curvatureMean = dataSum$curvature / (nRuns+1) 
    dataSum$ris = dataSum$ris / (nRuns+1) 
    dataSum$orientation = dataSum$orientation / (nRuns+1) 
  
    ind = which(dataSum$u < lT)
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
    dataSum$orientation  = dataSum$orientation[ind]
   
    ind = which(dataSum$v[ind] > uT)
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
    dataSum$orientation  = dataSum$orientation[ind]


    if(rm.na == T){
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
        dataSum$orientation  = dataSum$orientation[ind]

    }
    dataSum$curvature = sign(dataSum$cross)*acos(dataSum$angles)
    data[[i]] = dataSum

  }
  data
}




#get rim data for file i by circle fitting
get.inner.rim <- function(file){ 
  library(circular)
  x <- read.table( file, sep=",", header=F, colClasses="numeric")
  
  x <- x[1:(nrow(x)-1), ]
  circ = lsfit.circle(x)
  rim = list( mean = circ$coefficients[2:3], r=circ$coefficients[1], x=x)

   rim

}



get.outer.rim <- function(file){
  library(circular)

  x <- read.table( file, sep=",", header=F, colClasses="numeric")
  
  x <- x[1:(nrow(x)-1), ]
  circ = lsfit.circle(x)
  rim = list( mean = circ$coefficients[2:3], r=circ$coefficients[1], x=x)
  
  rim

}




#extract segments of length k from all features (read from the data with above
# methods)
extract.all.segments <- function(features, k, rm.jumps=k>1, offset=1 ){
  Slist = list()
  
  for(i in 1:length(features)){
    Slist[[i]] = extract.segments(X=features[[i]], k=k, rm.jumps=rm.jumps, offset=offset) 
    print(i)
  }

  Slist
}



#extract segements of length k from a single feature set X 
extract.segments  <- function(X, k, rm.jumps=k>1, offset = 1){    
  n = length( X$lengths )
    lengths=c()
    curvature=c()
    curvatureMean=c()
    dirim=c()
    dorim=c()
    time=c()
    x=c()
    y=c()
    rx=c()
    ry=c()
    orientation = c()
    for(j in 1:k){
      s = seq(j, n-k+j, by=offset)
      lengths = cbind( lengths, X$lengths[s])
        curvature = cbind( curvature, X$curvature[s])
        curvatureMean = cbind( curvatureMean, X$curvatureMean[s])
        dirim = cbind( dirim, X$dirim[s])
        dorim = cbind( dorim, X$dorim[s])
        time = cbind( time, X$time[s])
        x= cbind( x,X$x[s,1])
        y = cbind( y,X$x[s,2])
        rx= cbind( rx,X$ris[s,1])
        ry = cbind( ry,X$ris[s,2])
        orientation = cbind( orientation, X$orientation[s])

    }

  if(rm.jumps){
    jump = time[,k] - time[,1]
    ind  = which(jump == ( k-1))
    lengths = lengths[ind, ]
    curvature = curvature[ind, ]  
    curvatureMean = curvatureMean[ind, ]  
    dirim = dirim[ind, ]
    dorim = dorim[ind, ]
    time = time[ind, ]
    x = x[ind, ]
    y = y[ind, ]
    rx = rx[ind, ]
    ry = ry[ind, ]
    orientation = orientation[ind, ]
  }


  C = matrix(0, nrow=nrow(lengths), ncol=11)

    start1 = which(X$time > 5400 & X$dirim > 0)
    start1 = X$time[start1[1]]
    start2 = start1+300
    start3 = start2+600
    for( k in 1:nrow(C) ){
      C[k, 1] = sum( time[k, ] < start1 )
        C[k, 2] = sum( time[k, ] > start1 & time[k, ]  < start2 )
        C[k, 3] = sum( time[k, ] > start2 & time[k, ]  < start3 )
        C[k, 4] = sum( time[k, ] > start3 & time[k, ]  < 10800 )
        C[k, 5] = sum( time[k, ] > 10800 )
        C[k, 6] = which.max(C[k, 1:5])
        C[k, 7] = mean(time[k, ]) - start1

        C[k, 8] =  sum(dorim[k,] < 20)
        C[k, 9] = sum( dorim[k,] >= 20 & dirim[k, ] < -5)
        C[k, 10] = sum( dirim[k,] > -5 )
        C[k, 11] = which.max(C[k, 8:10])
    }
  colnames(C) <- c("before.odor", "odor1", "odor2", "odor3", "after.odor",
      "max.odor", "time.to.odor", "outer", "middler", "inner", "max.pos");

  list( lengths = lengths, curvature=curvature, curvatureMean=curvatureMean, dirim=dirim, dorim=dorim,
      time=time,x=x, y=y, rx=rx, ry=ry, C=as.data.frame(C), id=X$id,
      orientation=orientation)

}




#extract feature vectors from segments
extract.segment.features <- function(S){
  F = c()
  C = c()
  for(i in 1:nrow(S$lengths) ){
     
    ml = mean( S$lengths[i, ] )
    vl = var( S$lengths[i, ] )
    minl = min( S$lengths[i, ] )
    maxl = max( S$lengths[i, ] )

    mac = mean( abs( S$curvatureMean[i, ] ) )
    vc = var( S$curvatureMean[i, ] )
    minc = min( S$curvatureMean[i, ] )
    maxc = max( S$curvatureMean[i, ] )
    
    mo = mean( S$orientation[i,] )
    mao = mean( abs( S$orientation[i,])  )
    vo = var( S$orientation[i,] )
    mino = min( S$orientation[i, ] )
    maxo = max( S$orientation[i, ] )

    C = rbind(C, S$C[i, ])
    F = rbind(F, c(ml, vl, minl, maxl, mac, vc, minc, maxc, mo, vo, mino, maxo) )
  }

  colnames(F) <- c("mean.step", "var.step", "min.step", "max.step",
      "mean.abs.curvature", "var.curvature", "min.curvature", "max.curvature",
      "mean.orientation", "var.orientation", "min.orientation",
      "max.orientation" )

  list(F=F, C=C)
}


#extract feature vectors from segments
extract.all.segment.features <- function(Slist){
   Flist = list()
  
  for(i in 1:length(Slist)){
    Flist[[i]] = extract.segment.features(Slist[[i]]) 
    print(i)
  }

  Flist
}





#given the data X with n rows and conditions C (i.e. the matrix C from
#extract.segments) for each row extract submatrices condition on odor events
extract.condition.odor <- function(X, C){
  if( is.matrix(X) && ncol(X) > 1){
  Xbefore = X[C$max.odor == 1, ]
  Xduring1 = X[C$max.odor == 2 , ]
  Xduring2 = X[C$max.odor == 3 , ]
  Xduring3 = X[C$max.odor == 4 , ]
  Xafter = X[C$max.odor == 5, ]
  }
  else{
    Xbefore = X[ C$max.odor == 1]
    Xduring1 = X[ C$max.odor == 2 ]
    Xduring2 = X[ C$max.odor == 3 ]
    Xduring3 = X[ C$max.odor == 4 ]
    Xafter = X[ C$max.odor == 5 ]

  }
  list( Xbefore = Xbefore, 
        Xduring1 = Xduring1, Xduring2 = Xduring2,
        Xduring3 = Xduring3, Xafter = Xafter )

}

extract.condition.odor2 <- function(X, C){
  if( is.matrix(X) && ncol(X) > 1){
  Xbefore = X[C$max.odor == 1, ]
  Xduring = X[C$max.odor == 2 | C$max.odor == 3 | C$max.odor == 4 , ]
  Xafter = X[C$max.odor == 5, ]
  }
  else{
    Xbefore = X[ C$max.odor == 1]
    Xduring = X[C$max.odor == 2 | C$max.odor == 3 | C$max.odor == 4 ]
    Xafter = X[ C$max.odor == 5 ]

  }
  list( Xbefore = Xbefore, 
        Xduring = Xduring, 
        Xafter = Xafter )

}


extract.condition.odor.all <- function(Xlist, Slist){


  ans = NULL

  for(i in 1:length(Xlist) ){
    print(i)
    res = extract.condition.odor( Xlist[[i]], Slist[[i]]$C )
    if( is.null(ans) ){
      ans = res
    }
    else{
      for(i in 1:length(res)){
        ans[[i]] = rbind(ans[[i]], res[[i]])
      }
    }
  }

  ans
}









#Extract based on inside outside and odor period
extract.condition.odor.in.out <- function(X, C){

 
  BeforeInner = X[C$max.odor == 1 & C$max.pos == 3, ]
  BeforeMiddle = X[C$max.odor == 1 & C$max.pos == 2, ]
  BeforeOuter = X[C$max.odor == 1 & C$max.pos == 1, ]
  DuringInner = X[ (C$max.odor == 2 | C$max.odor == 3 | C$max.odor == 4 ) & C$max.pos == 3, ]
  DuringMiddle = X[ (C$max.odor == 2 | C$max.odor == 3 | C$max.odor == 4 ) & C$max.pos == 2, ]
  DuringOuter = X[ (C$max.odor == 2 | C$max.odor == 3 | C$max.odor == 4 ) & C$max.pos == 1, ]
  AfterInner = X[C$max.odor == 5 & C$max.pos == 3, ]
  AfterMiddle = X[C$max.odor == 5 & C$max.pos == 2, ]
  AfterOuter = X[C$max.odor == 5 & C$max.pos == 1, ]
  
  res = list(BeforInner = BeforeInner, BeforeMiddle = BeforeMiddle, BeforeOuter
      = BeforeOuter, DuringInner = DuringInner, DuringMiddle = DuringMiddle,
      DuringOuter = DuringOuter, AfterInner = AfterInner, AfterMiddle =
      AfterMiddle, AfterOuter =AfterOuter)
  res
}





extract.condition.odor.in.out.all <- function(Xlist, Slist){


  ans =NULL

  for(i in 1:length(Xlist) ){
    print(i)
    res = extract.condition.odor.in.out( Xlist[[i]], Slist[[i]]$C )
    if( is.null(ans) ){
      ans = res
    }
    else{
      for(i in 1:length(res)){
        ans[[i]] = rbind(ans[[i]], res[[i]])
      }
    }
  }

  ans
}




