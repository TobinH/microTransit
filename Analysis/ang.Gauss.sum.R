# experimenting with Tyler's (1987) maximum likelihood estimator for lambda
# of an angular central gaussian distribution for a 3d sphere
ang.Gauss.summ<-function(X, CI){
  lambda<-diag(1, nrow=3, ncol=3) # identity lamba matrix to start
  lambdahat<-matrix(data=0, nrow=3, ncol=3) # lambda estimate for Tyler's eqn (2)
  while (lambda[1,1] != lambdahat[1,1]){ # iterative procedure for ML estimate of lambda matrix
    sumli<-matrix(data=0, nrow=ncol(X), ncol=(ncol(X)))
    sumone<-vector(mode="numeric", length=1)
    sumone<-0
    for (i in 1:nrow(X)){
      sumli<-((X[i,]%*%t(X[i,]))/as.numeric(t(X[i,])%*%solve(lambda)%*%X[i,]))+sumli
      sumone<-(1/as.numeric(t(X[i,])%*%solve(lambda)%*%X[i,]))+sumone
    }
    lambda<-ncol(X)*sumli/sumone
    lambdahat<-ncol(X)/nrow(X)*sumli
  }
  gam<-eigen(lambda) # distribution axis weights and vectors
  
  aCI12<-asin(sqrt((nrow(X)*ncol(X))^-1
                         *(ncol(X)+2)^-1
                         *qchisq(1-(CI/2), df=ncol(X)-1)
                         *gam$values[2]
                         *gam$values[1])
                    /abs(gam$values[2]-gam$values[1])) # confidence interval from Axis I eigenvectors in Axis I-II plane, in radians
  aCI13<-asin(sqrt((nrow(X)*ncol(X))^-1
                         *(ncol(X)+2)^-1
                         *qchisq(1-(CI/2), df=ncol(X)-1)
                         *gam$values[3]
                         *gam$values[1])
                    /abs(gam$values[3]-gam$values[1])) # confidence interval from Axis I eigenvectors in Axis I-III plane, in radians
  aCI23<-asin(sqrt((nrow(X)*ncol(X))^-1
                         *(ncol(X)+2)^-1
                         *qchisq(1-(CI/2), df=ncol(X)-1)
                         *gam$values[3]
                         *gam$values[2])
                    /abs(gam$values[3]-gam$values[2])) # confidence interval around Axis I eigenvectors in Axis II-II plane, in radians
  aCI<-matrix(data=0, nrow=2, ncol=3)
  aCI[1,]<-c(aCI12, aCI13, aCI23)
  aCI[2,]<-aCI[1,]/pi*180
  res<-vector(mode="list", length=4)
  names(res)<-c("lambda.matrix", "lambda.eigval", "lambda.eigvec", "aCI")
  res$lambda.matrix<-lambda
  res$lambda.eigval<-gam$values
  names(res$lambda.eigval)<-c("Axis I", "Axis II", "Axis III")
  res$lambda.eigvec<-gam$vectors
  rownames(res$lambda.eigvec)<-colnames(X)
  colnames(res$lambda.eigvec)<-c("Axis I", "Axis II", "Axis III")
  res$aCI<-aCI
  colnames(res$aCI)<-c("I-II", "I-III", "II-III")
  rownames(res$aCI)<-c("radians", "degrees")
  return(res)
}
  
