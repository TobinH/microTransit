# binning tabulated data into generalized cortical dimensions for interspecific
# spatial comparison (for ISPH 2015). Follow this function with compliance matrix
# modeling.

# begin with 3 degree radial bins:

ISPH2015prep<-function(qPLMtab){
  results<-matrix(data=0,nrow=nrow(qPLMtab$pixels),ncol=10)
  # holds output
  results[,1:5]<-qPLMtab[,1:5]
  # transfer theta, phi, and xyz notation for pixel orientations
  qplMtab[,6:7]\<-scale(qPLMtab[,6:7],scale=FALSE,center=TRUE)
  # center pixel u,v positions
  results[,6]<-arctan(qPLMtab[,7]/qPLMtab[,6])%/%(2*pi*10)
  # place pixels in azimuth bins in 0.1 radian increments from the centroid
  results[,7]<-sqrt(qPLMtab[,7]^2+qPLMtab[,6]^2)
  #  get pixel distances from the centroid
  for (i in 0:120){
    dex<-which(results[,6]=i)
    halfdif<-(max(results[dex,7])-min(results[dex,7]))/2
    results[dex,7]<-trunc((results[dex,7]-(min(results[dex,7])+halfdif))/halfdif*5)
  }
  # scale pixel distances from centroid in each bin to (-5, 5)
  # columns six and seven are now equivalent distances in an idealized cylindrical
  # projection of the cortex
  
  # next model bending and torsion couples at body mass @ 1g--amplitudes aren't particularly vital here,
  # as this will only model differences in elastic behavior
  
  invisible(results)
  return(results)
}

