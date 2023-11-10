#' @title Place qPLMarr Data in Tabular Format
#'
#' @description \code{qPLMTabulate} casts data from the [u,v,3] array format of
#'   a qPLMarr object into a qPLMtab matrix.
#'
#' @param qPLMarr A qPLMarr object.
#'
#' @param low.pass A low-pass filter value (as an 8-bit grayscale level, 1-256)
#'   to remove noisy in-plane orientation data from low-theta pixels. Default is
#'   five.
#'
#' @return A [u*v,9] qPLMtab matrix, with columns (1) theta, (2) phi, (3-5) the
#'   Euclidean projection of theta and phi, (6-7) the [u,v] address of the pixel
#'   as an integer, and (8-9) the [u,v] position of the pixel in millimeters.
#'
#' @examples
#' #oldwd<-getwd()
#' #setwd(system.file("extdata", package = "microTransit"))
#' #load("testqPLMarr.R")
#' testqPLMtab<-qPLMTabulate(testqPLMarr)
#' testqPLMtab[1:5,]
#' save(testqPLMtab, file = "testqPLMtab.R")
#' #setwd(oldwd)
#'
#' @family qPLM Analysis Functions
#'
#' @export
#'


qPLMTabulate<-function(qPLMarr,low.pass=5/256){
  mashButton<-"y"
  if (attr(qPLMarr, "class") != "qPLMarr" | length(dim(qPLMarr))<3){
    mashButton<-readline(prompt = "This doesn't look like a qPLMarr object. Continue (y/n)? >")
  }
  if (mashButton !="y"){
    break
  }
  results<-vector(mode="list")
  nz.pix<-which(qPLMarr[,,2]>(low.pass), arr.ind=TRUE)
  # low pass retardance filter to remove "empty" pixels from analysis
  tabulated.data<-matrix(data=0,nrow=nrow(nz.pix),ncol=9)
  # setup by-pixel data matrix

  colnames(tabulated.data)<-c("theta","phi",
                              "x","y","z",
                              "u","v",
                              "u_mm", "v_mm")
  tabulated.data[,1]<-as.vector(qPLMarr[,,2][nz.pix]*pi/2)
  # pixel theta
  tabulated.data[,2]<-as.vector(qPLMarr[,,3][nz.pix]*pi)
  # pixel phi
  tabulated.data[,3]<-sin(tabulated.data[,1])*cos(tabulated.data[,2])
  # pixel orientation x
  tabulated.data[,4]<-sin(tabulated.data[,1])*sin(tabulated.data[,2])
  # pixel orientation y
  tabulated.data[,5]<-cos(tabulated.data[,1])
  # pixel orientation z
  tabulated.data[,6]<-nz.pix[,1]
  # pixel u position in image
  tabulated.data[,7]<-nz.pix[,2]
  # pixel v position in image
  tabulated.data[,8]<-nz.pix[,1]*attr(qPLMarr,"pix")/1000
  # pixel x position in image, calibrated to mm
  tabulated.data[,9]<-nz.pix[,2]*attr(qPLMarr,"pix")/1000
  # pixel y position in image, calibrated to to mm

  attr(tabulated.data, "thickness_um")<-attr(qPLMarr, "thickness_um")
  attr(tabulated.data, "wavelength_nm")<-attr(qPLMarr, "wavelength_nm")
  attr(tabulated.data, "birefringence")<-attr(qPLMarr, "birefringence")
  attr(tabulated.data, "pixel.size_um")<-attr(qPLMarr, "pixel.size_um")
  attr(tabulated.data, "ccw.skew_deg")<-attr(qPLMarr, "ccw.skew_deg")
  attr(tabulated.data, "dtype")<-attr(qPLMarr, "dtype")
  attr(tabulated.data, "class")<-"qPLMtab"

  return(tabulated.data)
}
