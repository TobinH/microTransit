#' @title Jones Matrix Representation of a Retarder Element
#'
#' @description \code{jonesRet} creates a Jones matrix for a retarder element,
#'   either with \enumerate{ \item (a) specified retardance or \item (b)
#'   specified birefringence, orientation, thickness, and wavelength,} at
#'   specified slow axis orientation.
#'
#' @details The Jones calculus provides a means to represent polarization of light with matrix and vector elements (Fowles 1989).
#'
#' @param retardance Shortcut to produce a specified retardance. Default is
#'   null.
#'
#' @param biref Birefringence of the material to be modeled. Required if
#'   \code{retardance} is null.
#'
#' @param theta.deg Out-of-plane orientation of the slow axis of the material.
#'   Required if \code{retardance} is null.
#'
#' @param d.um Thickness of the element in micrometers. Required if
#'   \code{retardance} is null.
#'
#' @param lambda Wavelength for simulation. Required if \code{retardance} is
#'   null.
#'
#' @param slow.phi.deg In-plane orientation of the slow axis of the material.
#'
#' @return A Jones matrix.
#'
#' @references Fowles, G. R. (1989). Introduction to modern optics. Courier Corporation.
#'
#' @family qPLM Simulation Functions
#'
#' @export

#
jonesRet<-function(retardance=NULL, biref, theta.deg, d.um, lambda, slow.phi.deg) {
  if (is.null(retardance)) {
    theta<-theta.deg/360*2*pi
    unit.ret<-biref*(sin(theta)^2)
    retardance<-unit.ret*d.um*1000/lambda
  }
  phi<-slow.phi.deg/360*2*pi
  out<-matrix(data=c(cos(retardance/2)+1i*sin(retardance/2)*cos(2*phi),
                     1i*sin(retardance/2)*sin(2*phi),
                     1i*sin(retardance/2)*sin(2*phi),
                     cos(retardance/2)-1i*sin(retardance/2)*cos(2*phi)),
              nrow=2, ncol=2)
  return(out)
}
