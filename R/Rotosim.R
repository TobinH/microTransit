#' @title Single-Pixel Simulation of Rotating Polarizer qPLM
#'
#' @description \code{Rotosim} uses Jones calculus to simulate the interactions
#' of a right circular polarizer, a birefringent specimen, and a rotating
#' analyzer in sequence, modeling the setup of a Glazer et al. (1996) method
#' qPLM scope.
#'
#' @param steps Integer; the number of analyzer positions to use. Default is
#'   six.
#'
#' @param sample a Jones matrix representation of a birefringent sample. Can be
#'   made using \code{jones.ret}.
#'
#' @return A list containing (1) |sin d|, a measure of retardance (varies with
#'   out-of-plane angle), (2-3) phi, the in-plane orientation of the specimen
#'   slow axis, and (4) Io, the transmissivity of the specimen.
#'
#' @references Glazer, A.M., Lewis, J.G., and Kaminsky, W., 1996. An automatic
#'   optical imaging system for birefringent media. \emph{Proc R Soc A
#'   452(1955)}: 2751 - 2765.
#'
#' @family qPLM Simulation Functions
#'
#' @export

# Single-pixel simulation of rotating polarizer qPLM methods
# Tobin Hieronymus 6-Oct-2014

Rotosim<-function(steps=6,sample){
  views<-matrix(data=NA, nrow=steps, ncol=4)
  RCP.light<-matrix(data=c(1/sqrt(2),(1/sqrt(2))*-1i), nrow=2, ncol=1)
  for (p in 1:steps){
    a.phi<-(p)*(pi/steps)
    analyzer<-matrix(data=c(cos(a.phi)^2,
                            sin(a.phi)*cos(a.phi),
                            cos(a.phi)*sin(a.phi),
                            sin(a.phi)^2), nrow=2, ncol=2)
    signal<-analyzer%*%sample%*%RCP.light
    views[p,1]<-Mod(signal[1,])^2+Mod(signal[2,])^2
    views[p,2]<-views[p,1]*sin(2*a.phi)
    views[p,3]<-views[p,1]*cos(2*a.phi)
    views[p,4]<-(a.phi/(2*pi))*360
  }
  a0<-sum(views[,1])/steps
  a1<-sum(views[,2])/(steps/2)
  a2<-sum(views[,3])/(steps/2)
  abs.sin.d<-(sqrt((a1^2)+(a2^2)))/a0
  phi1<-asin(a2/(sqrt(a1^2+a2^2)))/2
  phi1<-(phi1/(2*pi))*360
  phi2<-acos(-a1/(sqrt((a1^2)+(a2^2))))/2
  phi2<-(phi2/(2*pi))*360
  Io<-2*a0
  out<-list(abs.sin.d, phi1, phi2, Io)
}
