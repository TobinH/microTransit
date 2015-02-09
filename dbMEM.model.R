# dbMEM analysis function
# using spatial eigenfunctions from gen.sp.eigenfunction

dbMEM.model<-function(xyz,MEM){
  require(MASS)
  xnam <- paste0("x", 1:ncol(MEM))
  colnames(MEM)<-xnam
  xyznam<-c("x","y","z")
  colnames(xyz)<-xyznam
  fmlaup<-as.formula(paste(paste(xyznam,collapse= "+")," ~ ", paste(xnam, collapse= "+")))
  fmlalo<-as.formula(paste(paste(xyznam,collapse="+"),"~ 1"))                     
  fitupper<-lm(fmlaup)
  fitlower<-lm(fmlalo)
  results<-stepAIC(fit2,direction="forward",scope=list(upper=fit1,lower=fit2))
  return(results)
}
