# Read Rotopol (Werner Kaminsky) a1 bitmaps from known orientation materials and solve for birefringence.
# arguments:
# thickness: section thickness, in um
# wavelength: measurement wavelength, in nm
# angle: off-axis observation angle, _in degrees_ (0 is looking straight down the axis, 90 is flat to the slide)
# mask: Loads a black-and-white 8-bit bitmap that defines a region of interest. Black pixels (0) will be dropped, 
#   white pixels (255) will be kept. If mask=TRUE, a window will open to allow you to choose the mask file.

birefringence.estimate<-function(thickness, wavelength, angle, mask=FALSE) {
  require(bmp)
  require(EBImage)
  bmpfile<-choose.files(default="", caption="Select a1 bitmap", multi=FALSE)
  if (mask) {
    maskfile<-choose.files(default="", caption="Select mask bitmap")
  }
  sind<-read.bmp(bmpfile)
  sind<-sind[,,2]/255
  ifelse(mask, maskarray<-read.bmp(maskfile)/255, maskarray<-array(1,dim=dim(sind)))
  maskarray<-maskarray[,,2]
  sind<-sind*maskarray
  nz.pix<-which(sind>0, arr.ind=TRUE)
  result<-matrix(data=0,nrow=nrow(nz.pix),ncol=4)
  result[,1]<-(wavelength*asin(sind[nz.pix]))/(2*pi*thickness*1000*((sin(angle*2*pi/360))^2))
  result[,2]<-thickness
  result[,3]<-wavelength
  result[,4]<-angle
  colnames(result)<-c("birefringence", "thickness_um", "wavelength_nm", "angle_deg")
  return(result)
}