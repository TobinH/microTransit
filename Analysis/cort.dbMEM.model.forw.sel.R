# dbMEM analysis function--modeling with forward selection
# using spatial eigenfunctions from gen.sp.MEM
# arguments--
# qPLMtab: a tabulated qPLM dataset (output from qPLMtabulate.R)
# MEM: a Moran Eigenvector Map (output from gen.sp.MEM.R)
# penalty: IC penalty for adding parameters. 2 is AIC, log(n) is BIC

cort.dbMEM.model<-function(data,data.pos,MEM,penalty=2){
  require(qtlmt)
  pos.index<-match(data.frame(t(data.pos)), data.frame(t(MEM[,1:2])))
  # match eigenvectors to radial and centrifugal positions of pixel data
  eignf.data<<-as.data.frame(cbind(data,MEM[pos.index,3:ncol(MEM)]))
  # compose analysis data frame
  print(dim(eignf.data))
  print(colnames(eignf.data)[1:5])
  fmlaup<-as.formula(paste("cbind(sin_theta_doubled,cos_theta_doubled,sin_phi_corr,cos_phi_corr) ~ 1 + ", paste(colnames(MEM)[3:ncol(MEM)], collapse= "+")))
  fmlalo<-as.formula(cbind(sin_theta_doubled,cos_theta_doubled,sin_phi_corr,cos_phi_corr) ~ 1)                    
  fitupper<-lm(fmlaup,data=eignf.data)
  fitlower<-lm(fmlalo,data=eignf.data)
  results<-mStep(fitlower,direction="forward",scope=list(upper=fitupper,lower=fitlower),trace=TRUE,k=penalty)
  rm(eignf.data)
  return(results)
}
