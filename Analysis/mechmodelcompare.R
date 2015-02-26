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
  greatest.effect.b<-which.max(effect.comp[3,2,])
  highest.r.sq<-which.max(t(colSums(t(rsq.comp))))
  result<-list(coeff.comp=coeff.comp,rsq.comp=rsq.comp,effect.comp=effect.comp,bnum=greatest.effect.b,
               bmodel=link.mech.mod$modelfit[greatest.effect.b],
               rnum=highest.r.sq,rsqmodel=link.mech.mod$modelfit[highest.r.sq])
  return(result)
}

