# MEM modeling without variable selection--resulting in scalogram plot, fit, and anova
# arguments--
# MEM: a Moran Eigenvector Map matrix, output from gen.sp.MEM.R

dbMEM.model.scalogram<-function(MEM){
  eignf.data<-as.data.frame(MEM)
  fmla<-as.formula(paste("cbind(x,y,z) ~ 1 + ", paste(colnames(MEM)[5:ncol(MEM)], collapse= "+")))
  fit<-lm(fmla,data=eignf.data)
  anva<-anova(fit)
  results<-list(fit=NULL,anova=NULL)
  results$fit<-fit
  results$anova<-anva
  plot(results$anova[2:nrow(results$anova),3],type="l")
  return(results)
}
