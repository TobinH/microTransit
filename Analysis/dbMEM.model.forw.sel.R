# dbMEM analysis function--modeling with forward selection
# using spatial eigenfunctions from gen.sp.MEM
# arguments--
# qPLMtab: a tabulated qPLM dataset (output from qPLMtabulate.R)
# MEM: a Moran Eigenvector Map (output from gen.sp.MEM.R)
# penalty: IC penalty for adding parameters. 2 is AIC, log(n) is BIC

dbMEM.model<-function(qPLMtab,MEM,penalty=2){
  require(qtlmt)
  eignf.data<-as.data.frame(MEM)
  fmlaup<-as.formula(paste("cbind(x,y,z) ~ one + ", paste(xnam, collapse= "+")))
  fmlalo<-as.formula(cbind(x,y,z) ~ one)                    
  fitupper<-lm(fmlaup,data=eignf.data)
  fitlower<-lm(fmlalo,data=eignf.data)
  results<-mStep(fitlower,direction="forward",scope=list(upper=fitupper,lower=fitlower),trace=TRUE,k=penalty)
  return(results)
  }
}
