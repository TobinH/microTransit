gen.sp.dist<-function(x,cortical=FALSE) {
  require(ape)
  require(vegan)
  if(cortical) {
    posdist<-dist(x,method="manhattan")
    print(paste("Distance matrix completed at",Sys.time()))
    posdist[which(posdist>1)]<-4
  } else {
    posdist<-dist(x,method="euclidean")
    print(paste("Distance matrix completed at",Sys.time()))
    postree<-spantree(posdist)
    print(paste("Minimum spanning tree completed at",Sys.time()))
    posmax<-max(postree$dist)
    posdist[which(posdist>posmax)]<-posmax*4
  }
  centre <- function(D, n) {
    One <- matrix(1, n, n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  D<-as.matrix(posdist)
  n<-nrow(D)
  epsilon<-sqrt(.Machine$double.eps)
  delta1<-centre((-0.5*D^2),n)
  trace<-sum(diag(delta1))
  D.eig<-eigen(delta1)
  #min.eig<-min(D.eig$values)
  eig<-(D.eig$values)
  k<-length(which(eig>epsilon))
  #rel.eig<-eig[1:k]/trace
  #cum.eig<-cumsum(rel.eig)
  results<-sweep(D.eig$vectors[,1:k],2,sqrt(eig[1:k]),FUN="*")
  return(results)
}

#centre <- function(D, n) {
#  require(snow)
#  One <- matrix(1, n, n)
#  mat <- diag(n) - One/n
#  mat.cen <- parMM(cl,mat,(parMM(cl,D,mat)))
#}


#library(snow)
#library(parallel)
#library(doParallel)
#cl<-makeCluster(detectCores(),outfile="")
#registerDoParallel(cl)


#stopCluster(cl)
