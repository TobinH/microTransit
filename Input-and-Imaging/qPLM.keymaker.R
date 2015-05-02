# Rotopol.qPLM keymaker
# creates a key image for interpreting combined colatitude (theta) and longitude (phi)
# "_overview.tif" images

qPLM.keymaker<-function(pixel.size) {
  require(EBImage)
  require(colorspace)
  require(bmp)
  test.array<-array(data=NA, c(pixel.size,pixel.size,3))
  test.1<-matrix(0, nrow=length(test.array[,1,1]), ncol=length(test.array[1,,1]))
  x.center<-length(test.1[,1])/2
  y.center<-length(test.1[1,])/2
  for(i in 1:(dim(test.array)[1]/2)-1) {
    xx<-as.integer(x.center+i*cos(seq(0,2*pi,length.out=3600)))
    yy<-as.integer(y.center+i*sin(seq(0,2*pi,length.out=3600)))
    for(j in 1:3599) {
      test.1[[xx[j],yy[j]]]=(i/((dim(test.array)[1]/2)-1))
    }
  }
  test.2<-matrix(0, nrow=length(test.array[,1,1]), ncol=length(test.array[1,,1]))
  for(k in 1:nrow(test.2)){
    for(m in 1:ncol(test.2)){
      test.2[k,m]<-1-((atan((m-x.center)/(k-y.center))/(pi))+0.5)
    }
  }
  test.array[,,2]<-test.1
  test.array[,,3]<-test.2
  
  thetaseq<-seq(from=0, to=90)
  phiseq<-seq(from=0, to=180)
  pixLUT<-expand.grid(thetaseq, phiseq)
  rm(thetaseq, phiseq)
  # integer combinations of theta and phi for the look-up table
  
  PLUV.LUT<-polarLUV(pixLUT[,1]*73.2/90, pixLUT[,1]*57.65/90, pixLUT[,2]*360/180)
  # polar LUV colorspace encoding of each integer combination
  
  RGB.LUT<-as(PLUV.LUT,"RGB")
  rm(PLUV.LUT)
  # RGB value conversion of polar LUV look-up table
  
  pixmat<-as.matrix(pixLUT)
  rm(pixLUT)
  # n by 2 integer combination table
  
  RGBmat<-RGB.LUT@coords
  rm(RGB.LUT)
  gc()
  # n by 3 RGB look-up table
  
  immat<-cbind(as.vector(as.integer(test.array[,,2]*90)), as.vector((as.integer(test.array[,,3]*180))))
  # image pixel values m by 2 table
  
  LUindex<-match(data.frame(t(immat)), data.frame(t(pixmat)))
  rm(immat, pixmat)
  # look-up indices for image pixels
  
  key<-array(data=RGBmat[LUindex,], dim=dim(test.array))
  rm(RGBmat)
  rm(LUindex)
  gc()
  # cast looked-up RGB values into RGB image array
  
  output<-Image(key,colormode="Color")
  writeImage(output,file=paste(pixel.size,"_pixel_key.tif", sep=""), bits.per.sample=8L, type="tiff")
  # combined colatitude and azimuth color image out
  
  display(output)
  # pull up overview image in browser window
}