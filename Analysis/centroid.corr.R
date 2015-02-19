# setting orientation data to a centroid frame of reference for cortical bone
# resulting orientations are equivalent to 'unwrapping' the cortex

centroid.corr<-function(theta,phi,pos.x,pos.y) {
  positions<-cbind(pos.x,pos.y)
  positions<-scale(positions,center=TRUE,scale=FALSE)
  # centering pixel x,y positions to centroid
  corr.a<-cbind(theta,phi+atan(positions[,1]/positions[,2]))
  corr.o<-cbind(sin(corr.a[,1])*cos(corr.a[,2]),sin(corr.a[,1])*sin(corr.a[,2]),cos(corr.a[,1]))
  return(corr.o)
}
