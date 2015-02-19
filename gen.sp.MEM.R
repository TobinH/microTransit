# generate spatial eigenfunctions from qPLMtab data
# arguments--
# x: a qPLMtab$distance index (4x downsampled pixel position matrix)
# cortical: if TRUE, centers pixels and converts the cortex to a cylinder of standardized height and unit diameter
# returns an (n x p) matrix of p spatial eigenfunction scores for n pixels (Moran's Eigenvector Map)

gen.sp.MEM<-function(x,cortical=FALSE) {
  require(ape)
  require(vegan)
  pb<-winProgressBar(title=paste("gen.sp.MEM",format(Sys.time(), format="%H:%M")), min=0, max= 11, width=400)
  if(cortical) {
    setWinProgressBar(pb, 1, title=paste("gen.sp.MEM: centering cortex",format(Sys.time(), format="%H:%M")))
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
  setWinProgressBar(pb, 7, title=paste("gen.sp.MEM: eigenanalyzing delta1",format(Sys.time(), format="%H:%M")))
  D.eig<-eigen(delta1)
  #min.eig<-min(D.eig$values)
  eig<-(D.eig$values)
  setWinProgressBar(pb, 8, title=paste("gen.sp.MEM: retaining positive eigenfunctions",format(Sys.time(), format="%H:%M")))
  k<-length(which(eig>epsilon))
  #rel.eig<-eig[1:k]/trace
  #cum.eig<-cumsum(rel.eig)
  setWinProgressBar(pb, 9, title=paste("gen.sp.MEM: scaling eigenfunctions",format(Sys.time(), format="%H:%M")))
  results<-sweep(D.eig$vectors[,1:k],2,sqrt(eig[1:k]),FUN="*")
  setWinProgressBar(pb, 10, title=paste("gen.sp.MEM: removing NaN eignefunctions",format(Sys.time(), format="%H:%M")))
  if (any(is.nan(results))) {
    results<-results[,1:(which(is.nan(results),arr.ind=TRUE)[,2][1])-1]
  }
  setWinProgressBar(pb, 11, title=paste("gen.sp.MEM: done!",format(Sys.time(), format="%H:%M")))
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
