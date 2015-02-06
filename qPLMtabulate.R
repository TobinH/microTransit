# pull analysis data from microtransit object
# note that as of 4 Feb 15 microtransit is not a defined object class
# seems like a natural progression, but one thing at a time

qPLMtabulate<-function(x,low.pass=5/256){
  results<-list(pixels=NULL,distance=NULL)
  nz.pix<-which(x[,,2]>(low.pass), arr.ind=TRUE)
  # low pass retardance filter to remove "empty" pixels from analysis
  fourx.pos<-unique(nz.pix%/%4)
  # 4x4 binned pixel xy positions
  tabulated.data<-matrix(data=0,nrow=nrow(nz.pix),ncol=9)
  # setup by-pixel data matrix
  tabulated.data[,1]<-as.vector(x[,,2][nz.pix]*pi/2)
  # pixel theta
  tabulated.data[,2]<-as.vector(x[,,3][nz.pix]*pi)
  # pixel phi
  tabulated.data[,3]<-sin(tabulated.data[,1])*cos(tabulated.data[,2])
  # pixel orientation x
  tabulated.data[,4]<-sin(tabulated.data[,1])*sin(tabulated.data[,2])
  # pixel orientation y
  tabulated.data[,5]<-cos(tabulated.data[,1])
  # pixel orientation z
  tabulated.data[,6]<-nz.pix[,1]
  # pixel x position in image
  tabulated.data[,7]<-nz.pix[,2]
  # pixel y position in image
  tabulated.data[,8]<-nz.pix[,1]%/%4
  # pixel x position in 4x4 downsampled distances
  tabulated.data[,9]<-nz.pix[,2]%/%4
  # pixel y position in 4x4 downsampled distances
  results$pixels<-tabulated.data
  results$distance<-fourx.pos
  return(results)
}
