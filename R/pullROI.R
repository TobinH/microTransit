#' @title Interactive Sub-sampling of qPLM Data by Region of Interest
#'
#' @description \code{pullROI} sub-samples a \code{qPLMarr} object with an
#'   interactively selected rectangular ROI.
#'
#' @details no.
#'
#' @param qPLMarr A \code{qPLMarr} object.
#'
#' @param layer Value from 1-3; this layer will be displayed on-screen for
#'   selecting the ROI (1 = \strong{I}, 2 = \strong{theta}, 3 = \strong{phi}).
#'   The resulting ROI will include all of the 'stacked' layers from the
#'   \code{qPLMarr} object.
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
#' testROI<-pullROI(testqPLMarr, 2)
#'
#' @family qPLM Analysis Functions
#'
#' @export

pullROI<-function(qPLMarr, layer){
  #require(EBImage)
  mashButton<-"y"
  if (attr(qPLMarr, "class") != "qPLMarr"){
    mashButton<-readline(prompt = "This doesn't look like a qPLMarr object. Continue anyway (y/n)? > ")
    if (mashButton != "y"){
      break
    }
  }
  sel<-"n"
  while (sel != "y") {
    EBImage::display(qPLMarr[,,layer], method="raster")
    roi<-locator(2)
    rect(roi$x[1], roi$y[2], roi$x[2], roi$y[1], border="red")
    sub<-as.array(qPLMarr[roi$x[1]:roi$x[2],roi$y[1]:roi$y[2],])
    attributes(sub)<-c(attributes(sub), attributes(qPLMarr)[3:9])
    display(sub[,,layer], method="raster")
    sel<-readline("Keep ROI (y/n)? > ")
  }
  return(sub)
}
