# summary of xyz orientation data

scatter.matrix.summary<-function(qPLMtab){
  xyz<-as.matrix(qPLMtab$pixels[,6:7])
  results<-list(Tbar=NULL,R=NULL,eign=NULL,var=NULL)
  Tbar<-(t(xyz)%*%xyz)/length(xyz[,1])
  R<-sqrt(sum(colSums(xyz)^2))
  Tbar.Principal<-eigen(Tbar)
  results$Tbar<-Tbar
  # xyz var/covar matrix
  results$R<-R
  # resultant vector length
  results$eign<-Tbar.Principal
  # eigenvectors highlight major, middle, and minor axes in xyz space
  # eigenvalue comparison reveals information about distribution
  results$var<-1-(R/length(xyz[,1]))
  # Mardia's spherical variance
  return(results)
}
