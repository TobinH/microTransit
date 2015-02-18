# creates a jones matrix linear polarizer element at specified orientation phi
jones.pol<-function(phi) {
  rad<-phi/360*2*pi
  out<-matrix(data=c(cos(rad)^2, sin(rad)*cos(rad), cos(rad)*sin(rad), sin(rad)^2), nrow=2, ncol=2)
  return(out)
}
