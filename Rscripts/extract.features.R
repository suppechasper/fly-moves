
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


extract.all.features <- function(lT=6, rm.na=T){
  data <- list()
  for( i in 1:length(files)){
    print(i)
    data[[i]] = extract.features(i, lT, rm.na)
  }
  data
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




extract.all.segments <- function(features, k){
  Slist = list()
  
  for(i in 1:length(features)){
    Slist[[i]] = extract.segments(features[[i]], k) 
     print(i)
  }

  Slist
}



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
        C[k, 1] = sum( S$conds[[k]]$time < start1 )
          C[k, 2] = sum( time[k, ] > start1 & S$conds[[k]]$time < start2 )
          C[k, 3] = sum( time[k, ] > start2 & S$conds[[k]]$time < start3 )
          C[k, 4] = sum( time[k, ] > start3 & S$conds[[k]]$time < 10800 )
          C[k, 5] = sum( time[k,] > 10800 )
          C[k, 6] = which.max(C[k, 1:5])
          C[k, 7] = mean(time[k,]) - start1

          C[k, 8] =  sum(dorim[k,] < 20)
          C[k, 9] = sum( dorim[k,] >= 20 & S$conds[[k]]$dirim < -5)
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





