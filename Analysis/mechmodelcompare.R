# comparison of 'link function' mechanical models

mechmodelcompare<-function(link.mech.mod){
  coeff.comp<-array(data=0,dim=c(3,2,72))
  rsq.comp<-matrix(data=0,nrow=72,ncol=2)
  effect.comp<-array(data=0,dim=c(3,2,72))
  for (i in 1:72) {
    coeff.comp[,,i]<-link.mech.mod$modelfit[[i]]$coefficients
    rsq.comp[i,1]<-link.mech.mod$summ[[i]][[1]]$adj.r.squared
    rsq.comp[i,2]<-link.mech.mod$summ[[i]][[2]]$adj.r.squared
    effect.comp[,,i]<-link.mech.mod$modelfit[[i]]$effects[1:3,]
  }
  result<-list(coeff.comp=coeff.comp,rsq.comp=rsq.comp,effect.comp=effect.comp)
  return(result)
}

