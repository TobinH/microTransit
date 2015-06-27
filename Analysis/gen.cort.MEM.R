# generate spatial eigenfunctions for cortical comparison
# arguments--
# qPLMtab: a tabulated qPLM data set (output from qPLMtabulate)
# cortical: if TRUE, centers pixels and converts the cortex to a cylinder of standardized height and unit diameter
# maps: passes to options:nev in ARPACK, number of eigenvalues to compute
# returns an (n x p+4) matrix of (x,y,z) pixel orientation, a column of 1s, 
#   p spatial eigenfunction scores for n pixels (qPLM + Moran's Eigenvector Map)

gen.cort.MEM<-function(maps=300) {
  require(ape)
  require(vegan)
  require(igraph)
  cort<-expand.grid(1:360,-15:15)
  pb<-winProgressBar(title=paste("gen.cort.MEM", 
                                 format(Sys.time(), format="%H:%M")), 
                     min=0, max= 10, width=400)
  setWinProgressBar(pb, 1, title=paste("gen.cort.MEM: computing distances"
                                       ,format(Sys.time(), format="%H:%M")))
  posdist<-dist(cort,method="euclidean")
  setWinProgressBar(pb, 2, title=paste("gen.cort.MEM: minimum spanning tree"
                                       ,format(Sys.time(), format="%H:%M")))
  postree<-spantree(posdist)
  posmax<-max(postree$dist)
  posdist[which(posdist>posmax)]<-posmax*4
  centre <- function(D, n) {
    One <- matrix(1, n, n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  setWinProgressBar(pb, 3, title=paste("gen.cort.MEM: starting PCoA",
                                       format(Sys.time(), format="%H:%M")))
  D<-as.matrix(posdist)
  n<-nrow(D)
  epsilon<-sqrt(.Machine$double.eps)
  setWinProgressBar(pb, 4, title=paste("gen.cort.MEM: centering distance matrix"
                                       ,format(Sys.time(), format="%H:%M")))
  delta1<-centre((-0.5*D^2),n)
  setWinProgressBar(pb, 5, title=paste("gen.cort.MEM: trace of delta1"
                                       ,format(Sys.time(), format="%H:%M")))
  trace<-sum(diag(delta1))
  #setWinProgressBar(pb, 7, title=paste("gen.cort.MEM: eigenanalyzing delta1",format(Sys.time(), format="%H:%M")))
  #D.eig<-eigen(delta1)
  setWinProgressBar(pb, 6, 
                    title=paste("gen.cort.MEM: eigenanalyzing delta1 with ARPACK"
                                ,format(Sys.time(), format="%H:%M")))
  f2<-function(x, extra=NULL) { cat("."); as.vector(delta1 %*% x) }
  D.eig<-arpack(f2, sym=TRUE, 
                options=list(n=nrow(delta1),nev=maps,ncv=maps+2,
                             which="LM",maxiter=maps*1000))
  #min.eig<-min(D.eig$values)
  eig<-(D.eig$values)
  setWinProgressBar(pb, 7, 
                    title=paste("gen.cort.MEM: retaining positive eigenfunctions"
                                ,format(Sys.time(), format="%H:%M")))
  k<-length(which(eig>epsilon))
  #rel.eig<-eig[1:k]/trace
  #cum.eig<-cumsum(rel.eig)
  setWinProgressBar(pb, 8, title=paste("gen.sp.MEM: scaling eigenfunctions"
                                       ,format(Sys.time(), format="%H:%M")))
  evecs<-sweep(D.eig$vectors[,1:k],2,sqrt(eig[1:k]),FUN="*")
  setWinProgressBar(pb, 9, title=paste("gen.sp.MEM: removing NaN eigenfunctions"
                                       ,format(Sys.time(), format="%H:%M")))
  if (any(is.nan(evecs))) {
    evecs<-evecs[,1:(which(is.nan(evecs),arr.ind=TRUE)[,2][1])-1]
  }
  setWinProgressBar(pb, 10, 
                    title=paste("gen.cort.MEM: matching eigenfunctions to bins"
                                ,format(Sys.time(), format="%H:%M")))
  fnam <- paste("f_", 1:ncol(evecs),sep="")
  colnames(evec)<-fnam
  posnam<-c("radial","centrifugal")
  colnames(cort)<-posnam
  results<-as.matrix(cbind(cort,evecs))
  close(pb)
  print(paste("Finished gen.cort.MEM at",format(Sys.time(), format="%H:%M")))
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
