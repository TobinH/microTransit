# from qPLMtab to mechanical model estimates:

gen.mech.model<-function(qPLMtab){
  wgt.pts<-matrix(data=0,nrow=nrow(qPLMtab$pixels),ncol=3)
  pos<-matrix(data=0,nrow=nrow(qPLMtab$pixels),ncol=2)
  pos<-qPLMtab$pixels[,6:7]
  pos<-scale(pos,center=TRUE,scale=FALSE)
  posmod<-matrix(data=0,nrow=nrow(pos),ncol=3)
  posmod[,1:2]<-pos[,1:2]
  posmod[,3]<-0
  wgt.pts[,1]<-sin(2*qPLMtab$pixels[,1])*cos(qPLMtab$pixels[,2])*
    (sqrt(pos[,1]^2+pos[,2]^2))
  # x in cartesian spherical projection of CFO weighted by distance from 
  # centroid--theta out to pi from 1/2pi
  wgt.pts[,2]<-sin(2*qPLMtab$pixels[,1])*sin(qPLMtab$pixels[,2])*
    (sqrt(pos[,1]^2+pos[,2]^2))
  # y in ""
  wgt.pts[,3]<-cos(2*qPLMtab$pixels[,1])*(sqrt(pos[,1]^2+pos[,2]^2))
  # z in ""
  #Tw<-(t(wgt.pts)%*%wgt.pts)/length(wgt.pts[,1])
  results<-list(A=NULL,mod.points=NULL)
  # results$eigTw<-eigen(Tw)
  Smod<-t(posmod)%*%wgt.pts/nrow(qPLMtab$pixels)
  dec.Smod<-svd(Smod)
  A<-dec.Smod$u%*%t(dec.Smod$v)
  results$A<-A
  points<-t(A%*%t(posmod))
  points<-points[,]/sqrt(points[,1]^2+points[,2]^2+points[,3]^2)
  results$mod.points<-points
  return(results)
  # troubleshooting results for now
  # form of results should be a per-pixel xyz estimate
}
  
