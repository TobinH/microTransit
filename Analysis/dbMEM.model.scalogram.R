# MEM modeling without variable selection--resulting in scalogram plot, fit, and anova
# arguments--
# qPLMtab: a tabulated qPLM dataset (output from qPLMtabulate.R)
# MEM: a Moran Eigenvector Map matrix, output from gen.sp.MEM.R

dbMEM.model.scalogram<-function(qPLMtab,MEM){
  xnam <- paste("x_", 1:ncol(MEM),sep="")
  xyznam<-c("x","y","z")
  xyz<-matrix(data=qPLMtab$pixels[,3:5],ncol=3)
  colnames(xyz)<-xyznam
  one<-matrix(data=1,nrow=nrow(xyz),ncol=1)
  colnames(one)<-"one"
  pixindex<-match(as.data.frame(t(qPLMtab$pixels[,8:9])),as.data.frame(t(qPLMtab$distance)))
  MEMpix<-matrix(data=0,nrow=nrow(xyz),ncol=ncol(MEM))
  MEMpix<-MEM[pixindex,]
  colnames(MEMpix)<-xnam
  eignf.data<<-as.data.frame(cbind(xyz,one,MEMpix))
  fmla<-as.formula(paste("cbind(x,y,z) ~ one + ", paste(xnam, collapse= "+")))
  fit<-lm(fmla,data=eignf.data)
  anva<-anova(fit)
  results<-list(fit=NULL,anova=NULL)
  results$fit<-fit
  results$anova<-anva
  plot(anva[2:nrow(anva),3],type="l")
  return(results)
}
