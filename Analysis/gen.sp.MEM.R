# generate spatial eigenfunctions from qPLMtab data
# arguments--
# qPLMtab: a tabulated qPLM data set (output from qPLMtabulate)
# cortical: if TRUE, centers pixels and converts the cortex to a cylinder of standardized height and unit diameter
# maps: passes to options:nev in ARPACK, number of eigenvalues to compute
# returns an (n x p+4) matrix of (x,y,z) pixel orientation, a column of 1s, 
#   p spatial eigenfunction scores for n pixels (qPLM + Moran's Eigenvector Map)

gen.sp.MEM<-function(qPLMtab,cortical=FALSE,maps=300) {
  require(ape)
  require(vegan)
  require(igraph)
  pb<-winProgressBar(title=paste("gen.sp.MEM",format(Sys.time(), format="%H:%M")), min=0, max= 12, width=400)
  if(cortical) {
    setWinProgressBar(pb, 1, title=paste("gen.sp.MEM: centering cortex",format(Sys.time(), format="%H:%M")))
    x<-qPLMtab$distance
    x<-scale(x,center=TRUE,scale=TRUE)
    rad<-matrix(data=0,nrow=nrow(x),ncol=3)
    rad[,1]<-sqrt(x[,1]^2+x[,2]^2)
    rad[,2]<-x[,1]/rad[,1]
    rad[,3]<-x[,2]/rad[,1]
    rad[,1]<-scale(rad[,1],center=TRUE,scale=TRUE)
    setWinProgressBar(pb, 2, title=paste("gen.sp.MEM: computing distances",format(Sys.time(), format="%H:%M")))
    posdist<-dist(rad,method="euclidean")
    setWinProgressBar(pb, 3, title=paste("gen.sp.MEM: minimum spanning tree",format(Sys.time(), format="%H:%M")))
    postree<-spantree(posdist)
    posmax<-max(postree$dist)
    posdist[which(posdist>posmax)]<-posmax*4
    } else {
    setWinProgressBar(pb, 2, title=paste("gen.sp.MEM: computing distances",format(Sys.time(), format="%H:%M")))
    posdist<-dist(x,method="euclidean")
    setWinProgressBar(pb, 3, title=paste("gen.sp.MEM: minimum spanning tree",format(Sys.time(), format="%H:%M")))
    postree<-spantree(posdist)
    posmax<-max(postree$dist)
    posdist[which(posdist>posmax)]<-posmax*4
  }
  centre <- function(D, n) {
    One <- matrix(1, n, n)
    mat <- diag(n) - One/n
    mat.cen <- mat %*% D %*% mat
  }
  setWinProgressBar(pb, 4, title=paste("gen.sp.MEM: starting PCoA",format(Sys.time(), format="%H:%M")))
  D<-as.matrix(posdist)
  n<-nrow(D)
  epsilon<-sqrt(.Machine$double.eps)
  setWinProgressBar(pb, 5, title=paste("gen.sp.MEM: centering distance matrix",format(Sys.time(), format="%H:%M")))
  delta1<-centre((-0.5*D^2),n)
  setWinProgressBar(pb, 6, title=paste("gen.sp.MEM: trace of delta1",format(Sys.time(), format="%H:%M")))
  trace<-sum(diag(delta1))
  #setWinProgressBar(pb, 7, title=paste("gen.sp.MEM: eigenanalyzing delta1",format(Sys.time(), format="%H:%M")))
  #D.eig<-eigen(delta1)
  setWinProgressBar(pb, 7, title=paste("gen.sp.MEM: eigenanalyzing delta1 with ARPACK",format(Sys.time(), format="%H:%M")))
  f2<-function(x, extra=NULL) { cat("."); as.vector(delta1 %*% x) }
  D.eig<-arpack(f2, sym=TRUE, options=list(n=nrow(delta1),nev=maps,ncv=maps+2,which="LM",maxiter=maps*1000))
  #min.eig<-min(D.eig$values)
  eig<-(D.eig$values)
  setWinProgressBar(pb, 8, title=paste("gen.sp.MEM: retaining positive eigenfunctions",format(Sys.time(), format="%H:%M")))
  k<-length(which(eig>epsilon))
  #rel.eig<-eig[1:k]/trace
  #cum.eig<-cumsum(rel.eig)
  setWinProgressBar(pb, 9, title=paste("gen.sp.MEM: scaling eigenfunctions",format(Sys.time(), format="%H:%M")))
  evecs<-sweep(D.eig$vectors[,1:k],2,sqrt(eig[1:k]),FUN="*")
  setWinProgressBar(pb, 10, title=paste("gen.sp.MEM: removing NaN eignefunctions",format(Sys.time(), format="%H:%M")))
  if (any(is.nan(evecs))) {
    evecs<-evecs[,1:(which(is.nan(evecs),arr.ind=TRUE)[,2][1])-1]
  }
  setWinProgressBar(pb, 11, title=paste("gen.sp.MEM: matching eigenfunctions to pixels",format(Sys.time(), format="%H:%M")))
  xnam <- paste("x_", 1:ncol(evecs),sep="")
  xyznam<-c("x","y","z")
  xyz<-matrix(data=qPLMtab$pixels[,3:5],ncol=3)
  colnames(xyz)<-xyznam
  one<-matrix(data=1,nrow=nrow(xyz),ncol=1)
  colnames(one)<-"one"
  pixindex<-match(as.data.frame(t(qPLMtab$pixels[,8:9])),as.data.frame(t(qPLMtab$distance)))
  MEMpix<-matrix(data=0,nrow=nrow(xyz),ncol=ncol(evecs))
  MEMpix<-evecs[pixindex,]
  colnames(MEMpix)<-xnam
  results<-as.matrix(cbind(xyz,one,MEMpix))
  setWinProgressBar(pb, 12, title=paste("gen.sp.MEM: done!",format(Sys.time(), format="%H:%M")))
  close(pb)
  print(paste("Finished gen.sp.MEM at",format(Sys.time(), format="%H:%M")))
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
