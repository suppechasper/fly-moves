
##
#frame odor on:   5400
#frame odor off: 10801
#fly dimensions 25 x 12 pixels
#instantaneous velocity < 3 pixels /s , values abouve 4 p/s probably errorenous
#peak fluctation is 0.15 pixels, sd 0.09 pixels


files = c(
"120202video12_xypts_transformed.csv",
"120202video18_xypts_transformed.csv",
"120202video6_xypts_transformed.csv",
"120209video6_xypts_transformed.csv",
"120218video13_xypts_transformed.csv",
"120218video7_xypts_transformed.csv",
"120316video5_xypts_transformed.csv",
"120323video6_xypts_transformed.csv",
"120328video7_xypts_transformed.csv",
"120418video7_xypts_transformed.csv",
"120427video7_xypts_transformed.csv",
"120503video6_xypts_transformed.csv",
"120517video16_xypts_transformed.csv",
"120530video8_xypts_transformed.csv",
"120605video38_xypts_transformed.csv",
"120607video6_xypts_transformed.csv",
"120612video36_xypts_transformed.csv",
"120613video6_xypts_transformed.csv",
"120614video8_xypts_transformed.csv",
"120808video13_xypts_transformed.csv",
"120809video7_xypts_transformed.csv",
"120828video6_xypts_transformed.csv",
"120829video1_xypts_transformed.csv",
"120911video8_xypts_transformed.csv",
"120913video7_xypts_transformed.csv",
"120920video8_xypts_transformed.csv",
"120925video6_xypts_transformed.csv",
"120928video6_xypts_transformed.csv",
"130213video8_xypts_transformed.csv",
"130214video9_xypts_transformed.csv"
)




rimFiles <- c(
"120202transformedrimpoints.matinner_transformed.csv",
"120202transformedrimpoints.matinner_transformed.csv",
"120202transformedrimpoints.matinner_transformed.csv",
"120209transformedrimpoints.matinner_transformed.csv",
"120218transformedrimpoints.matinner_transformed.csv",
"120218transformedrimpoints.matinner_transformed.csv",
"120316transformedrimpoints.matinner_transformed.csv",
"120323transformedrimpoints.matinner_transformed.csv",
"120328transformedrimpoints.matinner_transformed.csv",
"120418transformedrimpoints.matinner_transformed.csv",
"120427transformedrimpoints.matinner_transformed.csv",
"120503transformedrimpoints.matinner_transformed.csv",
"120517transformedrimpoints.matinner_transformed.csv",
"120530transformedrimpoints.matinner_transformed.csv",
"120605transformedrimpoints.matinner_transformed.csv",
"120607transformedrimpoints.matinner_transformed.csv",
"120612transformedrimpoints.matinner_transformed.csv",
"120613transformedrimpoints.matinner_transformed.csv",
"120614transformedrimpoints.matinner_transformed.csv",
"120808transformedrimpoints.matinner_transformed.csv",
"120809transformedrimpoints.matinner_transformed.csv",
"120828transformedrimpoints.matinner_transformed.csv",
"120829transformedrimpoints.matinner_transformed.csv",
"120911transformedrimpoints.matinner_transformed.csv",
"120913transformedrimpoints.matinner_transformed.csv",
"120920transformedrimpoints.matinner_transformed.csv",
"120925transformedrimpoints.matinner_transformed.csv",
"120928transformedrimpoints.matinner_transformed.csv",
"130213transformedrimpoints.matinner_transformed.csv",
"130214transformedrimpoints.matinner_transformed.csv"
)



outerRimFiles <- c(
"120202transformedrimpoints.matouter_transformed.csv",
"120202transformedrimpoints.matouter_transformed.csv",
"120202transformedrimpoints.matouter_transformed.csv",
"120209transformedrimpoints.matouter_transformed.csv",
"120218transformedrimpoints.matouter_transformed.csv",
"120218transformedrimpoints.matouter_transformed.csv",
"120316transformedrimpoints.matouter_transformed.csv",
"120323transformedrimpoints.matouter_transformed.csv",
"120328transformedrimpoints.matouter_transformed.csv",
"120418transformedrimpoints.matouter_transformed.csv",
"120427transformedrimpoints.matouter_transformed.csv",
"120503transformedrimpoints.matouter_transformed.csv",
"120517transformedrimpoints.matouter_transformed.csv",
"120530transformedrimpoints.matouter_transformed.csv",
"120605transformedrimpoints.matouter_transformed.csv",
"120607transformedrimpoints.matouter_transformed.csv",
"120612transformedrimpoints.matouter_transformed.csv",
"120613transformedrimpoints.matouter_transformed.csv",
"120614transformedrimpoints.matouter_transformed.csv",
"120808transformedrimpoints.matouter_transformed.csv",
"120809transformedrimpoints.matouter_transformed.csv",
"120828transformedrimpoints.matouter_transformed.csv",
"120829transformedrimpoints.matouter_transformed.csv",
"120911transformedrimpoints.matouter_transformed.csv",
"120913transformedrimpoints.matouter_transformed.csv",
"120920transformedrimpoints.matouter_transformed.csv",
"120925transformedrimpoints.matouter_transformed.csv",
"120928transformedrimpoints.matouter_transformed.csv",
"130213transformedrimpoints.matouter_transformed.csv",
"130214transformedrimpoints.matouter_transformed.csv"
)




#read file number index and extract features
extract.features <- function(index, lT=6, rm.na=T, std=0){

  x <- read.table(files[index], sep=",", header=T)


  x = x + matrix(rnorm(nrow(x)*ncol(x), sd=std), ncol = ncol(x))

  irim = get.inner.rim(index);
  orim = get.outer.rim(index);
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



#Extract features with noise added and avergae over nRuns
extract.all.features.expected <- function(lT=6, rm.na=T, nRuns=10, std=0.06){
  data <- list()
  for( i in 1:length(files)){
    print(i)
    dataSum = extract.features(i, lT, rm.na=F)
    dataSum$curvatureOrig = dataSum$curvature
    for(k in 1:nRuns){
      dataTmp = extract.features(i, lT, rm.na=F, std)
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
get.inner.rim <- function(i){
  library(circular)
  x <- read.table(rimFiles[i], sep=",", header=F)
  
  x <- x[1:(nrow(x)-1), ]
  circ = lsfit.circle(x)
  rim = list( mean = circ$coefficients[2:3], r=circ$coefficients[1], x=x)

   rim

}



get.outer.rim <- function(i){
  library(circular)

  x <- read.table(outerRimFiles[i], sep=",", header=F)
  
  x <- x[1:(nrow(x)-1), ]
  circ = lsfit.circle(x)
  rim = list( mean = circ$coefficients[2:3], r=circ$coefficients[1], x=x)
  
  rim

}







#--- older methods

extract.all.features <- function(lT=6, rm.na=T){
  data <- list()
  for( i in 1:length(files)){
    print(i)
    data[[i]] = extract.features(i, lT, rm.na)
  }
  data
}





extract.segments <- function(index, len, off = len/2){
  x <- read.table(files[index], sep=",", header=T)
  
  irim = get.inner.rim(index);
  orim = get.outer.rim(index);
 
  v1 = x;
  v1[,1] = v1[,1] - irim$mean[1] 
  v1[,2] = v1[,2] - irim$mean[2] 
  dirim = irim$r - sqrt( rowSums(v1^2) )


  v2 = x;
  v2[,1] = v2[,1] - orim$mean[1] 
  v2[,2] = v2[,2] - orim$mean[2] 
  dorim = orim$r - sqrt( rowSums(v2^2) )

  
    
  d = dim(x)
  segsX1 = matrix(0, nrow=d[1]-len, ncol=len)
  segsX2 = matrix(0, nrow=d[1]-len, ncol=len)
  segsDO = matrix(0, nrow=d[1]-len, ncol=len)
  segsDI = matrix(0, nrow=d[1]-len, ncol=len)
  times = matrix(0, nrow=d[1]-len, ncol=len)
  for( i in seq(1,(d[1]-len), off) ){
    segsX1[i, ] = x[i:(i+len-1), 1]    
    segsX2[i, ] = x[i:(i+len-1), 2] 
    segsDI[i, ] = dirim[i:(i+len-1)]   
    segsDO[i, ] = dorim[i:(i+len-1)]
    times[i, ] = i:(i+len-1)   
  }
  list(segsX1 = segsX1, segsX2 = segsX2, times=times, dorim=segsDO, dirim =
      segsDI)
}




extract.conditions <- function(data){
  conds <- c()
  for( i in 1:length(files)){
    print(i)
    conds = rbind(conds, cbind( data[[i]]$dirim, data[[i]]$dorim, data[[i]]$time, i) )
  }
  colnames(conds) <-  c("dirim", "dorim", "time", "id") 

  as.data.frame(conds)

}




extract.all.segments <- function(len, off = len/2){
  data <- list()
  for( i in 1:length(files)){
    print(i)
    data[[i]] = extract.segments(i, len, off);
  }
  data
}



clean.segments <- function(data){
   for( i in 1:length(data)){
    print(i)
    x = data[[i]]
    n <- nrow(x$segsX1)
    dx1 = (x$segsX1[1:(n-1), ] - x$segsX1[2:n, ])^2
    dx2 = (x$segsX2[1:(n-1), ] - x$segsX2[2:n, ])^2
    d = apply(dx1+dx2, 1, max)
    ind = which( d < 40 )
    x$segsX1 = x$segsX1[ind, ]
    x$segsX2 = x$segsX2[ind, ]
    x$dirim  = x$dirim[ind, ]
    x$dorim  = x$dorim[ind, ]
    x$times  = x$times[ind, ]
    data[[i]] = x
  }
 
  data
}



extract.segments.conditions <- function( data ){
  conds <- c()
  for( i in 1:length(files)){
    print(i)
    conds = rbind(conds, cbind( rowMeans(data[[i]]$dirim),
          rowMeans(data[[i]]$dorim), rowMeans( data[[i]]$times ), i) )
  }
  colnames(conds) <-  c("dirim", "dorim", "time", "id") 

  as.data.frame(conds)
}




extract.segment.features.sorted <- function(data, l, off=l/3){

  index = seq(1, length(data$lengths)-l-1, off)
  n =length(index)
  fl = matrix( 0, nrow=n, ncol=l)
  fa = matrix( 0, nrow=n, ncol=l)
  fc = matrix( 0, nrow=n, ncol=l)
  
  fdi = matrix( 0, nrow=n, ncol=1)
  fdo = matrix( 0, nrow=n, ncol=1)
  ft = matrix( 0, nrow=n, ncol=1)

  for(k in 1:n ){
    i = index[k]
    fl[k,] = sort( data$lengths[ i:(i+l-1) ] , na.last=T)
    fa[k,] = sort( data$angles[ i:(i+l-1) ], na.last=T)
    fc[k,] = sort( data$cross[ i:(i+l-1) ], na.last=T)

    fdi[k,1] = mean( data$dirim[ i:(i+l-1) ], na.rm=T)
    fdo[k,1] = mean( data$dorim[ i:(i+l-1) ], na.rm=T)
    ft[k,1] = mean( data$time[ i:(i+l-1) ], na.rm=T)
  }

  features <- list(lengths=fl, angles=fa, cross = fc,seg.length=l, time=ft,
      id=data$id, dorim=fdo, dirim = fdi)

}


extract.all.segment.features.sorted <- function(data, l, off = l/3){

  segs <- c()
  for( i in 1:length(data)){
    print(i)
    seg =  extract.segment.features.sorted( data[[i]], l, off)
    segs  = rbind(segs, cbind(seg$angles, seg$cross, seg$lengths, seg$time, 
          seg$dorim, seg$dirim, seg$id ) ) 
  }
  colnames(segs) = c(rep("angles", l), rep("cross", l), rep("lengths", l),
      "time", "dorim", "dirim", "id" ) 
  segs

}


extract.segment.features <- function(data, l, off=l/3){

  index = seq(1, length(data$lengths)-l-1, off)
  n =length(index)
  fl = matrix( 0, nrow=n, ncol=l)
  fa = matrix( 0, nrow=n, ncol=l)
  fc = matrix( 0, nrow=n, ncol=l)
  
  fdi = matrix( 0, nrow=n, ncol=1)
  fdo = matrix( 0, nrow=n, ncol=1)
  ft = matrix( 0, nrow=n, ncol=1)

  for(k in 1:n ){
    i = index[k]
    fl[k,] =  data$lengths[ i:(i+l-1) ]
    fa[k,] =  data$angles[ i:(i+l-1) ]
    fc[k,] =  data$cross[ i:(i+l-1) ]

    fdi[k,1] = mean( data$dirim[ i:(i+l-1) ], na.rm=T)
    fdo[k,1] = mean( data$dorim[ i:(i+l-1) ], na.rm=T)
    ft[k,1] = mean( data$time[ i:(i+l-1) ], na.rm=T)
  }

  features <- list(lengths=fl, angles=fa, cross = fc,seg.length=l, time=ft,
      id=data$id, dorim=fdo, dirim = fdi)

}


extract.all.segment.features <- function(data, l, off = l/3){

  segs <- c()
  for( i in 1:length(data)){
    print(i)
    seg =  extract.segment.features( data[[i]], l, off)
    segs  = rbind(segs, cbind(seg$angles, seg$cross, seg$lengths, seg$time, 
          seg$dorim, seg$dirim, seg$id ) ) 
  }
  colnames(segs) = c(rep("angles", l), rep("cross", l), rep("lengths", l),
      "time", "dorim", "dirim", "id" ) 
  segs

}





extract.segment.features.summary <- function(data, l){

  fl = matrix( 0, nrow=length(data$lengths)-l+1, 3)
  fa = matrix( 0, nrow=length(data$lengths)-l+1, 3)
  fc = matrix( 0, nrow=length(data$lengths)-l+1, 3)
  
  fdi = matrix( 0, nrow=length(data$lengths)-l+1, 1)
  fdo = matrix( 0, nrow=length(data$lengths)-l+1, 1)
  ft = matrix( 0, nrow=length(data$lengths)-l+1, 1)

  for(i in 1:nrow(fa) ){
    fl[i,1] = mean( data$lengths[ i:(i+l-1) ] , na.rm=T)
    fl[i,2] = median( data$lengths[ i:(i+l-1) ] , na.rm=T)
    fl[i,3] = stats::sd( data$lengths[ i:(i+l-1) ] , na.rm=T)
    fl[i,4] = max( data$lengths[ i:(i+l-1) ] , na.rm=T)

    fa[i,1] = mean( data$angles[ i:(i+l-1) ] , na.rm=T)
    fa[i,2] = median( data$angles[ i:(i+l-1) ] , na.rm=T)
    fa[i,3] = stats::sd( data$angles[ i:(i+l-1) ] , na.rm=T)

    fc[i,1] = mean( data$cross[ i:(i+l-1) ], na.rm=T)
    fc[i,2] = median( data$cross[ i:(i+l-1) ], na.rm=T)
    fc[i,3] = stats::sd( data$cross[ i:(i+l-1) ], na.rm=T )

    fdi[i,1] = mean( data$dirim[ i:(i+l-1) ], na.rm=T)
    fdo[i,1] = mean( data$dorim[ i:(i+l-1) ], na.rm=T)
    ft[i,1] = mean( data$time[ i:(i+l-1) ], na.rm=T)
  }

  features <- list(lengths=fl, angles=fa, cross = fc,seg.length=l, time=ft,
      id=data$id, dorim=fdo, dirim = fdi)

}


extract.all.segment.features.summary <- function(data, l){

  segs <- c()
  for( i in 1:length(data)){
    print(i)
    seg =  extract.segment.features( data[[i]], l)
    segs  = rbind(segs, cbind(seg$angles, seg$cross, seg$lengths, seg$time, 
          seg$dorim, seg$dirim, seg$id ) ) 
  }
  colnames(segs) = c(rep("angles", 3), rep("cross", 3), rep("lengths", 4),
      "time", "dorim", "dirim", "id" ) 
  segs

}


