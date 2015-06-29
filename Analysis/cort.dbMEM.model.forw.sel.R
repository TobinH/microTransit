# dbMEM analysis function--modeling with forward selection
# using spatial eigenfunctions from gen.sp.MEM
# arguments--
# qPLMtab: a tabulated qPLM dataset (output from qPLMtabulate.R)
# MEM: a Moran Eigenvector Map (output from gen.sp.MEM.R)
# penalty: IC penalty for adding parameters. 2 is AIC, log(n) is BIC

cort.dbMEM.model<-function(data,data.pos,MEM,penalty=2){
  require(qtlmt)
  pos.index<-match(data.frame(t(data.pos)), data.frame(t(MEM[,1:2])))
  
  eignf.data<-as.data.frame(cbind(data,MEM[pos.index,]))
  colnames(eignf.data[,1:4])<-c("st","ct","sp","cp")
  colnames(eignf.data[,5:ncol(MEM)+2])<-colnames(MEM[,3:ncol(MEM)])
  fmlaup<-as.formula(paste("cbind(st,ct,sp,cp) ~ one + ", paste(fnam, collapse= "+")))
  fmlalo<-as.formula(cbind(st,ct,sp,cp) ~ one)                    
  fitupper<-lm(fmlaup,data=eignf.data)
  fitlower<-lm(fmlalo,data=eignf.data)
  results<-mStep(fitlower,direction="forward",scope=list(upper=fitupper,lower=fitlower),trace=TRUE,k=penalty)
  return(results)
}
