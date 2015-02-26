# alternative 'link function' linear model for bending + torsional stress and collagen fiber orientation
# using Identity function I() to explicitly cast the model terms
# multiple iterations of the model around the range of possible neutral axes
# AIC to select 'best' model

link.mech.model<-function(qPLMtab) {
  theta<-qPLMtab$pixels[,1]
  phi<-qPLMtab$pixels[,2]
  v<-qPLMtab$pixels[,6]
  w<-qPLMtab$pixels[,7]
  tresist<-abs(cos(phi-atan2(w,v)))*cos(theta)
  bresist<-sin(2*theta)
  mechdata<-as.data.frame(cbind(tresist,bresist,v,w))
  mechdata[,3:4]<-scale(mechdata[,3:4],scale=FALSE,center=TRUE)
  modelfit<-vector("list", length=72)
  modelANOVA<-vector("list", length=72)
  summ<-vector("list", length=72)
  for (i in seq(from=0, to=355, by=5)) {
    a<-i/180*pi
    mechfmla<-as.formula(cbind(tresist,bresist) ~ 1 + I(sqrt(v^2+w^2)) + I(cos(a)*((v*tan(a))-w)))
    fit<-lm(mechfmla,data=mechdata)
    anva<-anova(fit)
    modelfit[[i/5+1]]<-fit
    modelANOVA[[i/5+1]]<-anva
    summ[[i/5+1]]<-summary(fit)
  }
  result<-list(modelfit=modelfit,modelANOVA=modelANOVA,summ=summ)
  #result$best.fit<-modelwheel[which.min(modelAIC)]
  return(result)
}
