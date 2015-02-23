# MEM modeling without variable selection--resulting in scalogram plot, fit, and anova
# arguments--
# MEM: a Moran Eigenvector Map matrix, output from gen.sp.MEM.R

dbMEM.model.scalogram<-function(MEM){
  eignf.data<-as.data.frame(MEM)
  fmla<-as.formula(paste("cbind(x,y,z) ~ one + ", paste(colnames(MEM)[4:ncol(MEM)], collapse= "+")))
  fit<-lm(fmla,data=eignf.data)
  anva<-anova(fit)
  results<-list(fit=NULL,anova=NULL)
  results$fit<-fit
  results$anova<-anva
  #plot(anva[2:nrow(anva),3],type="l", par(mar=c(0,0,0,0)))
  #par(mar=c(5,4,4,2)+0.1)
  return(results)
}
