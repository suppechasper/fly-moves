save.data <- function(Slist, fly.type="WT"){
 X <- c()

 for(i in 1:length(Slist) ){
   
   S <- Slist[[i]]
   x <-  cbind(S$x, S$y, S$rx, S$ry, S$lengths, S$curvature, S$curvatureMean,
       S$orientation, S$dirim, S$dorim, S$time, S$C, rep(S$id, length(S$x)) )
 
   colnames(x) <- c( "x", "y", "rx", "ry", "step.length", "curvature",
     "curvatureMean", "orientation", "distance.inner.rim", "distance.outer.rim",
     "time", colnames(Slist[[1]]$C), "fly.id" )

   X <- rbind(X, x)
 
 }

#colnames(X) <- c( "x", "y", "rx", "ry", "step.length", "curvature",
#     "curvatureMean", "orientation", "distance.inner.rim", "distance.outer.rim",
#     colnames(Slist[[1]]$C), "fly.id" )

  df = as.data.frame(X)
  df$fly.type = fly.type

  df
}



