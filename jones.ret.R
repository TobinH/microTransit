# creates a jones matrix for a retarder element, either with specified retardance
# or with specified birefringence, orientation, thickness, and wavelength,
# at specified slow axis orientation
jones.ret<-function(retardance=NULL, biref, theta.deg, d.um, lambda, slow.phi.deg) {
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
