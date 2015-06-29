# individual taxon-element matrix prep for ISPH full analysis
# notes for remaining steps:
#      1. concatenate matrices by element with a taxon flag
#      2.feed into pVARPART 

ISPH2015dataprep<-function(taxon.short,
                           element.abbr,
                           N,S,W,E,
                           cw=TRUE,
                           mass,
                           element.J,
                           pixel=7.5832259) {
  qPLMobj<-Rotopol.qPLM.script(paste(taxon.short, element.abbr, sep="."),
                               paste("*", element.abbr, "*full*", sep=""),
                               paste("*", element.abbr, "*mask*", sep=""),
                               N, S, W, E,
                               pics.out=FALSE,
                               pixel=pixel)
  qPLMt<-qPLMtabulate(qPLMobj, low.pass=1/256)
  qPLMr<-resolve.theta.ambiguity(qPLMt, cw=cw)
  qPLMe<-elastic.model(qPLMr, mass, element.J)
  qPLMs<-spatial.model(qPLMt,120,10)
  angles<-matrix(data=0,nrow=nrow(qPLMt$pixels),ncol=4)
  angles[,1]<-sin(2*qPLMt$pixels[,1])
  angles[,2]<-cos(2*qPLMt$pixels[,1])
  angles[,3]<-sin(qPLMr[,2]-atan2(qPLMe[,2], qPLMe[,1]))
  angles[,4]<-cos(qPLMr[,2]-atan2(qPLMe[,2], qPLMe[,1]))
  result<-cbind(angles,qPLMe,qPLMs[,6:7])
  colnames(result)<-c("sin_theta*2", "cos_theta*2",
                      "sin_phi_corr","cos_phi_corr",
                      "u","v",
                      "iso.raw","ortho.raw","diff.raw",
                      "iso.scaled","ortho.scaled","diff.scaled",
                      "spatial.radial","spatial.centrifugal")
  invisible(result)
  return(result)
}