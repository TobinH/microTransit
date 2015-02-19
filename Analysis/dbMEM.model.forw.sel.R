# dbMEM analysis function--modeling with forward selection
# using spatial eigenfunctions from gen.sp.MEM
# arguments--
# qPLMtab: a tabulated qPLM dataset (output from qPLMtabulate.R)
# MEM: a Moran Eigenvector Map (output from gen.sp.MEM.R)
# regress: if FALSE, this function matches pixels to their MEM score without running a model
# penalty: IC penalty for adding parameters. 2 is AIC, log(n) is BIC

dbMEM.model<-function(qPLMtab,MEM,regress=TRUE,penalty=2){
  require(qtlmt)
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
  if (!regress) {
    return(eignf.data[,5:ncol(eignf.data)])
  } else {
  fmlaup<-as.formula(paste("cbind(x,y,z) ~ one + ", paste(xnam, collapse= "+")))
  fmlalo<-as.formula(cbind(x,y,z) ~ one)                    
  fitupper<-lm(fmlaup,data=eignf.data)
  fitlower<-lm(fmlalo,data=eignf.data)
  results<-mStep(fitlower,direction="forward",scope=list(upper=fitupper,lower=fitlower),trace=TRUE,k=penalty)
  return(results)
  }
}
