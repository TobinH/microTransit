# binning tabulated data into generalized cortical dimensions for interspecific
# spatial comparison (for ISPH 2015). Follow this function with compliance matrix
# modeling.


spatial.model<-function(qPLMtab,rad.bin,half.centr.bin){
  results<-matrix(data=0,nrow=nrow(qPLMtab$pixels),ncol=10)
  # holds output
  results[,1:5]<-qPLMtab$pixels[,1:5]
  # transfer theta, phi, and xyz notation for pixel orientations
  qPLMtab$pixels[,6:7]<-scale(qPLMtab$pixels[,6:7],scale=FALSE,center=TRUE)
  # center pixel u,v positions
  results[,6]<-ceiling((atan2(qPLMtab$pixels[,7],qPLMtab$pixels[,6])+pi)
                       /(2*pi)*rad.bin)
  # place pixels in radial (azimuth) bins from the centroid
  results[,7]<-sqrt((qPLMtab$pixels[,7]^2)+(qPLMtab$pixels[,6]^2))
  #  get pixel distances from the centroid
  for (i in 1:rad.bin){
    dex<-which(results[,6]==i)
    halfdif<-(max(results[dex,7])-min(results[dex,7]))/2
    results[dex,7]<-trunc((results[dex,7]-(min(results[dex,7])+halfdif))
                          /halfdif*half.centr.bin)
  }
  # scale pixel distances from centroid 
  # in each bin to (-half.centr.bin, half.centr.bin)
  # columns six and seven are now equivalent posiions in an idealized cylindrical
  # projection of the cortex with n radial bins and 2p+1 centrifugal bins
  
  invisible(results)
  return(results)
}

