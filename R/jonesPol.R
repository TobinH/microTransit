#' @title Jones Matrix Representation of a Linear Polarizer Element
#'
#' @description \code{jonesPol} creates a Jones matrix linear polarizer element
#' at specified orientation phi.
#'
#' @details none yet.
#'
#' @param phi in-plane angle for the axis of the polarizer element.
#'
#' @return A Jones matrix.
#'
#' @family qPLM Simulation Functions
#'
#' @export

jonesPol<-function(phi) {
  rad<-phi/360*2*pi
  out<-matrix(data=c(cos(rad)^2, sin(rad)*cos(rad), cos(rad)*sin(rad), sin(rad)^2), nrow=2, ncol=2)
  return(out)
}
