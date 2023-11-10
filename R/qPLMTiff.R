#' @title Produce Overview TIFFs From qPLM Arrays
#'
#' @description \code{qPLMTiff} creates greyscale .tiff images of out-of-plane
#'   orientation (theta), in-plane orientation (phi), transmittance (trans), and
#'   a false-color representation of combined orientation data (overview).
#'
#' @details The false-color image of the "overview" \code{.tiff} uses the
#'   LUV colorspace to produce a perceptually similar brightness for equal
#'   levels of theta, and encodes phi with different hues.
#'
#' @param sample.name An identifier for the sample, used as a stem for the
#'   output file names.
#'
#' @param qPLMarr A qPLMarr object, e.g., output from \code{buildqPLM}.
#'
#' @param theta.scale Boolean; if TRUE, automatically scales "_overview.tif"
#'   image color scale to peak at maximum measured theta value (default is
#'   FALSE).
#'
#' @param theta.max Numerical; if theta.scale is FALSE, specifies an upper end
#'   (degrees) for "_overview.tif" color scale (default is 90 degrees).
#'
#' @param invert.theta.scale Boolean; if TRUE, bright colors point toward
#'   viewer, brightness falls off as orientations go to parallel with the slide.
#'   Useful for sections with predominantly out-of-plane fibers (default is
#'   FALSE).
#'
#' @return Null; .tiff files saved to working directory.
#'
#' @family qPLM Illustration Functions
#'
#' @export

qPLMTiff<-function(sample.name,
                   qPLMarr,
                   theta.scale=FALSE,
                   theta.max=90,
                   invert.theta.scale=FALSE){
  mashButton<-"y"
  if (attr(qPLMarr, "class") != "qPLMarr"){
    mashButton<-readline(prompt = "This doesn't look like a qPLMarr object. Continue anyway (y/n)? > ")
    if (mashButton != "y"){
      break
    }
  }
  if(theta.scale){
    theta.max<-ceiling(max(qPLMarr[,,2]*90))
  }
  theta.max<-ceiling(theta.max)
  if (invert.theta.scale) {
    thetaseq<-seq(from=90, to=90-theta.max)
  }
  else{
    thetaseq<-seq(from=0, to=theta.max)
  }
  phiseq<-seq(from=0, to=180)
  pixLUT<-expand.grid(thetaseq, phiseq)
  rm(thetaseq, phiseq)
  # integer combinations of theta and phi for the look-up table

  print("Building colorspace: LUV encoding")
  PLUV.LUT<-colorspace::polarLUV((pixLUT[,1]*0.75+22.5)*90/theta.max,
                                 pixLUT[,1]*57.65/theta.max,
                                 pixLUT[,2]*360/180)
  # polar LUV colorspace encoding of each integer combination

  print("Building colorspace: RGB translation")
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

  print("Color mapping")
  if (invert.theta.scale){
    theta.values<-as.vector(as.integer((1-qPLMarr[,,2])*90))
  }
  else{
    theta.values<-as.vector(as.integer(qPLMarr[,,2]*90))
  }

  immat<-cbind(theta.values, as.vector((as.integer(qPLMarr[,,3]*180))))
  # image pixel values m by 2 table

  LUindex<-match(data.frame(t(immat)), data.frame(t(pixmat)))
  rm(immat, pixmat)
  # look-up indices for image pixels

  coded<-array(data=RGBmat[LUindex,], dim=dim(qPLMarr))
  rm(RGBmat)
  rm(LUindex)
  gc()
  # cast looked-up RGB values into RGB image array

  print("Writing image: _trans")
  trans<-EBImage::Image(qPLMarr[,,1], colormode="Grayscale")
  EBImage::writeImage(trans, file=paste(sample.name,"_trans.tif", sep=""), bits.per.sample=8L, type="tiff")
  # Transmission (I) grayscale image out

  print("Writing image: _theta")
  theta<-Image(qPLMarr[,,2], colormode="Grayscale")
  writeImage(theta, file=paste(sample.name, "_theta.tif", sep=""), bits.per.sample=8L, type="tiff")
  # colatitude (theta) grayscale image out

  print("Writing image: _phi")
  phi<-Image(qPLMarr[,,3], colormode="Grayscale")
  EBImage::writeImage(phi, file=paste(sample.name, "_phi.tif",sep=""), bits.per.sample=8L, type="tiff")
  rm(trans, theta, phi)
  gc()
  # azimuth (phi) grayscale image out

  print("Writing image: _overview")
  output<-EBImage::Image(coded,colormode="Color")
  EBImage::writeImage(output,file=paste(sample.name,"_overview.tif", sep=""), bits.per.sample=8L, type="tiff")
  # combined colatitude and azimuth color image out

}
