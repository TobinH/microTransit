#' @title Interactive Sub-sampling of qPLM Data by Pixel count
#'
#' @description \code{pullqPLMSample} sub-samples a specified number of pixels
#'   from a \code{qPLMarr} object around an interactively selected point.
#'
#' @param qPLMarr A \code{qPLMarr} object.
#'
#' @param layer Value from 1-3; this layer will be displayed on-screen for
#'   selecting the ROI (1 = \strong{I}, 2 = \strong{theta}, 3 = \strong{phi}).
#'   The resulting subsample will include all of the 'stacked' layers from the
#'   \code{qPLMarr} object.
#'
#' @param pixels approximate number of pixels required for the sample (rounding
#'   for a square region may pull pixels more than specified by this argument).
#'
#' @return A \code{qPLMarr} object with the subsampled pixels. This object will
#'   maintain all the attributes of the original.
#'
#' @examples
#' #oldwd<-getwd()
#' #setwd(system.file("extdata", package = "microTransit"))
#' #load("testqPLMarr.R")
#' #setwd(oldwd)
#'
#' testSample<-pullqPLMSample(testqPLMarr)
#'
#' @family qPLM Analysis Functions
#'
#' @export


pullqPLMSample<-function(qPLMarr, layer, pixels){
  #require(EBImage)
  mashButton<-"y"
  if (attr(qPLMarr, "class") != "qPLMarr"){
    mashButton<-readline(prompt = "This doesn't look like a qPLMarr object. Continue anyway (y/n)? > ")
    if (mashButton != "y"){
      break
    }
  }
  frameSize<-ceiling(sqrt(pixels)/2)
  sel<-"n"
  while (sel != "y") {
    EBImage::display(qPLMarr[,,layer], method="raster")
    roi<-locator(1)
    rect((roi$x)-frameSize, (roi$y)-frameSize, (roi$x)+frameSize, (roi$y)+frameSize, border="red")
    sub<-as.array(qPLMarr[((roi$x)-frameSize):((roi$x)+frameSize), ((roi$y)-frameSize):((roi$y)+frameSize),])
    attributes(sub)<-c(attributes(sub), attributes(qPLMarr)[3:9])
    display(sub[,,layer], method="raster")
    sel<-readline("Keep subsample (y/n)? > ")
  }
  return(sub)
}
