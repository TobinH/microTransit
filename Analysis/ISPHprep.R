# binning tabulated data into generalized cortical dimensions for interspecific
# spatial comparison (for ISPH 2015). Follow this function with compliance matrix
# modeling.

# begin with 3 degree radial bins:

spatial.model<-function(qPLMtab){
  results<-matrix(data=0,nrow=nrow(qPLMtab$pixels),ncol=10)
  # holds output
  results[,1:5]<-qPLMtab$pixels[,1:5]
  # transfer theta, phi, and xyz notation for pixel orientations
  qPLMtab$pixels[,6:7]<-scale(qPLMtab$pixels[,6:7],scale=FALSE,center=TRUE)
  # center pixel u,v positions
  results[,6]<-ceiling((atan2(qPLMtab$pixels[,7],qPLMtab$pixels[,6])+pi)
                       /(2*pi)*360)
  # place pixels in azimuth bins in 1 degree increments from the centroid
  results[,7]<-sqrt((qPLMtab$pixels[,7]^2)+(qPLMtab$pixels[,6]^2))
  #  get pixel distances from the centroid
  for (i in 1:360){
    dex<-which(results[,6]==i)
    halfdif<-(max(results[dex,7])-min(results[dex,7]))/2
    results[dex,7]<-trunc((results[dex,7]-(min(results[dex,7])+halfdif))
                          /halfdif*15)
  }
  # scale pixel distances from centroid in each bin to (-15, 15)
  # columns six and seven are now equivalent posiions in an idealized cylindrical
  # projection of the cortex with 360 radial bins and 31 centrifugal bins
  
  invisible(results)
  return(results)
}

