#' @title Set Tangential Frame of Reference for Cortical Bone Axial Data
#'
#' @description \code{centroidCorr} changes the in-plane orientation data of a
#'   \code{qPLMtab} object from global reference frame to orientation relative
#'   to a tangent from the centroid of the object cross-section. Resulting
#'   orientations are in a sense 'unwrapping' the cortex around the centroid.
#'
#' @details For some questions (e.g., resistance of cortical bone to torsion),
#'   orientation in plane relative to a centroid is informative. This function
#'   returns orientation relative to a 'naive' assumption of the structural
#'   centroid of a cross-section.
#'
#' @param qPLMtab qPLMtab object to process.
#'
#' @return A \code{qPLMtab} object; a matrix with columns \enumerate{ \item (1)
#'   theta, \item (2) 'corrected' phi, \item (3-5) the 'unwrapped' Euclidean x,
#'   y, and z axis coordinates, and \item (6-7) the original pixel u,v
#'   positions.}
#'
#' @examples
#' #oldwd<-getwd()
#' #setwd(system.file("extdata", package = "microTransit"))
#' #load("testqPLMtab.R")
#' #setwd(oldwd)
#'
#' testCorrtab<-centroidCorr(testqPLMtab)
#'
#' @family qPLM Analysis Functions
#'
#' @export

centroidCorr<-function(qPLMtab) {
  mashButton<-"y"
  if (attr(qPLMtab, "class") != "qPLMtab"){
    mashButton<-readline(prompt = "This doesn't look like a qPLMtab object. Continue anyway (y/n)? > ")
    if (mashButton != "y"){
      break
    }
  }
  theta<-qPLMtab[,1]
  phi<-qPLMtab[,2]
  positions<-qPLMtab[,6:7]
  positions<-scale(positions,center=TRUE,scale=FALSE)
  # centering pixel x,y positions to centroid
  corr.a<-cbind(theta,phi-atan(positions[,2]/positions[,1])+pi/2)
  corr.o<-cbind(sin(corr.a[,1])*cos(corr.a[,2]),sin(corr.a[,1])*sin(corr.a[,2]),cos(corr.a[,1]))
  results<-cbind(corr.a, corr.o, qPLMtab[,6:7])
  attr(results, "thickness_um")<-attr(qPLMtab, "thickness_um")
  attr(results, "wavelength_nm")<-attr(qPLMtab, "wavelength_nm")
  attr(results, "birefringence")<-attr(qPLMtab, "birefringence")
  attr(results, "pixel.size_um")<-attr(qPLMtab, "pixel.size_um")
  attr(results, "ccw.skew_deg")<-attr(qPLMtab, "ccw.skew_deg")
  attr(results, "dtype")<-paste("centroid-corrected", attr(qPLMtab, "dtype"))
  attr(results, "class")<-"qPLMtab"
  return(results)
}
