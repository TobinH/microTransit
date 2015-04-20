# Alternative 'link function' linear model for bending + torsional 
# stress and collagen fiber orientation,  using Identity function 
# I() to explicitly cast the model terms.
# Multiple iterations of the model around the range of possible neutral axes.
# #AIC to select 'best' model3->not implemented
# global note: I should change v,w pixel notation to standard u,v
# pixel notation to avoid confusion

link.mech.model<-function(qPLMtab) {
  theta<-qPLMtab$pixels[,1]
  phi<-qPLMtab$pixels[,2]
  v<-qPLMtab$pixels[,6]
  w<-qPLMtab$pixels[,7]
  # load relevant variables from tabulated qPLM data object
  tresist<-abs(cos(phi-atan2(w,v)))*cos(theta)
  # simple model of torsional modulus around centroid
  bresist<-sin(2*theta)
  # simple model of compressive/tensile modulus
  mechdata<-as.data.frame(cbind(tresist,bresist,v,w))
  mechdata[,3:4]<-scale(mechdata[,3:4],scale=FALSE,center=TRUE)
  # combine modulus models and centered pixel v,w values
  modelfit<-vector("list", length=72)
  modelANOVA<-vector("list", length=72)
  summ<-vector("list", length=72)
  # objects to receive results of iterative fit process
  for (i in seq(from=0, to=355, by=5)) {
    a<-i/180*pi
    mechfmla<-as.formula(cbind(tresist,bresist) 
                         ~ 1 + I(sqrt(v^2+w^2)) + I(cos(a)*((v*tan(a))-w)))
    fit<-lm(mechfmla,data=mechdata)
    anva<-anova(fit)
    modelfit[[i/5+1]]<-fit
    modelANOVA[[i/5+1]]<-anva
    summ[[i/5+1]]<-summary(fit)
  }
  # fit model to data with 5-degree increments of neutral bending axis
  result<-list(modelfit=modelfit,modelANOVA=modelANOVA,summ=summ)
  # 
  #result$best.fit<-modelwheel[which.min(modelAIC)]
  return(result)
}
