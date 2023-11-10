#' @title Create an 3D Skin for Interpreting qPLM Overview Images
#'
#' @description \code{qPLMSphereSkinner} creates a skin image for placing
#' colatitude (theta) and longitude (phi) on a sphere (e.g., in Adobe
#' Illustrator).
#'
#' @details the \code{.tiff} image produced by this function can be used to skin
#'   a hemispherical object for inclusion in 3D-modeled figures.
#'
#' @param pixel.size Dimension in pixels for square output image.
#'
#' @param theta.max Maximum value of theta. Used when the qPLM object to
#'   interpret has a truncated theta range--take this value directly from the
#'   "theta.max" attribute of the qPLM object.
#'
#' @return Silent return in R. Saves an image with the name "x_pixel_skin.tif"
#'   where x is the value of pixel.size, and displays the image to a browser
#'   window.
#'
#' @family qPLM Illustration Functions
#'
#' @export


# Rotopol.qPLM sphere skinner
# creates a skin image for placing colatitude (theta) and longitude (phi) on a sphere

qPLMSphereSkinner<-function(pixel.size, theta.max) {
  #require(EBImage)
  #requireNamespace(colorspace)
  #require(bmp)

  colatitude.limit<-as.integer(pixel.size*theta.max/180)

  test.array<-array(data=NA, c(pixel.size, colatitude.limit,3))

  test.1<-matrix(0, nrow=pixel.size, ncol=colatitude.limit)
  for(i in 1:colatitude.limit){
    test.1[,i]<-i/(pixel.size/2)
  }

  test.2<-matrix(0, nrow=pixel.size, ncol=colatitude.limit)
  for(j in 1:pixel.size){
    test.2[j,]<-j/pixel.size
  }

  test.array[,,2]<-test.1
  test.array[,,3]<-test.2

  thetaseq<-seq(from=0, to=theta.max)
  phiseq<-seq(from=0, to=180)
  pixLUT<-expand.grid(thetaseq, phiseq)
  rm(thetaseq, phiseq)
  # integer combinations of theta and phi for the look-up table

  PLUV.LUT<-colorspace::polarLUV((pixLUT[,1]*0.75+22.5)*theta.max/90, pixLUT[,1]*57.65/theta.max, pixLUT[,2]*360/180)
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

  immat<-cbind(as.vector(as.integer(test.array[,,2]*theta.max)), as.vector((as.integer(test.array[,,3]*180))))
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
  EBImage::writeImage(output,file=paste(pixel.size, "_by_", theta.max, "_pixel_skin.tif", sep=""), bits.per.sample=8L, type="tiff")
  # combined colatitude and azimuth color image out

  EBImage::display(output)
  # pull up overview image in browser window
}
