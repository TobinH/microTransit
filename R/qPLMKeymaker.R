#' @title Create an Orientation Key for Interpreting qPLM Overview Images
#'
#' @description \code{qPLMKeymaker} creates a key image for interpreting
#' combined colatitude (theta) and longitude (phi) "_overview.tif" images.
#'
#' @details The \code{.tiff} file produced by this function should appear as a
#'   colorful circle on a black background. Imagine this image as a view of the
#'   northern hemisphere of a globe, with the prime meridian directed to the right.
#'   A ray projecting from the center of the globe through the north pole
#'   (colatitude 0°) would emerge in the center of the circle, and be color-coded
#'   black. A ray emerging with a colatitude of 45°and a longitude of 45°would
#'   be color-coded as a half-saturated blue. A ray emerging at the equator with
#'   a longitude of 135° would be color-coded as a fully saturated green.
#'
#' @param pixel.size Dimension in pixels for square output image.
#'
#' @return Silent return in R. Saves an image with the name "x_pixel_key.tif"
#'   where x is the value of pixel.size, and displays the image to a browser
#'   window.
#'
#' @family qPLM Illustration Functions
#'
#' @export
#'


qPLMKeymaker<-function(pixel.size) {
  #require(EBImage)
  #requireNamespace(colorspace)
  #require(bmp)
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

  PLUV.LUT<-colorspace::polarLUV((pixLUT[,1]*0.75+22.5), pixLUT[,1]*57.65/90, pixLUT[,2]*360/180)
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

  output<-EBImage::Image(key,colormode="Color")
  EBImage::writeImage(output,file=paste(pixel.size,"_pixel_key.tif", sep=""), bits.per.sample=8L, type="tiff")
  # combined colatitude and azimuth color image out

  EBImage::display(output)
  # pull up overview image in browser window
}
