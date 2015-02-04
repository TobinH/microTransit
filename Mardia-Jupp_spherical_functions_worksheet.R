# spherical stats for qPLM data from Mardia and Jupp 2000, Upton and Fingleton 1985.
# Tobin Hieronymus, 02 Feb 2015
# thieronymus@neomed.edu

# starting point--read.Rotopol.R function output, an i by j by 3 dimensional array, x and y are pixel positions.
# orientations are restricted axial--i.e., reflections of theta are identical, antipodal projections of theta and phi are identical
# [,,2] contains angular estimates of theta in the range (0,1), representing real space range of (0,pi/2) radians
# [,,3] contains angular measurements of phi in the range (0,1), representing real space range of (0,pi) radians

# initial questions: what's the distribution of orientation in a region of interest?
# what's the distribution of orientations with respect to the long axis of a long bone?

# first approach: pull array into a matrix of n non-zero-pixels by (x,y, and z of orientation data + pixel position i,j)

# using code from earlier circular stats functions
qPLMtabulate<-function(x,low.pass=5){
  nonzero.pixels<-which(x[,,2]>(low.pass), arr.ind=TRUE)
  angular.data<-cbind(x[cbind(nonzero.pixels,2)]/512*pi,x[cbind(nonzero.pixels,3)]/256*pi)
  raw.orientations<-cbind(sin(angular.data[,1])*cos(angular.data[,2]),sin(angular.data[,1])*sin(angular.data[,2]),cos(angular.data[,1]))
  positions<-scale(nonzero.pixels,center=TRUE,scale=FALSE)
  centr.corr.angular<-cbind(angular.data[,1],angular.data[,2]+atan(positions[,1]/positions[,2]))
  centr.corr.orientations<-cbind(sin(centr.corr.angular[,1])*cos(centr.corr.angular[,2]),sin(centr.corr.angular[,1])*sin(centr.corr.angular[,2]),cos(centr.corr.angular[,1]))
  Inertia<-t(positions)%*%positions
  Tbar.n<-(t(raw.orientations)%*%raw.orientations)
  Tbar<-Tbar.n/length(raw.orientations[,1])
  Tbar.corr.n<-t(centr.corr.orientations)%*%centr.corr.orientations
  Tbar.corr<-Tbar.corr.n/length(centr.corr.orientations[,1])
  R.raw<-sqrt(sum(colSums(raw.orientations)^2))
  R.corr<-sqrt(sum(colSums(centr.corr.orientations)^2))
  Xsxn.principal<-eigen(Inertia)
  Tbar.principal<-eigen(Tbar)
  Tbar.corr.principal<-eigen(Tbar.corr)
}

# notes for remaining steps 2-Feb-15:
# per pixel matrix of x,y,z points D, then scatter matrix T=D'D--eigenanalyze for principal axes
# per pixel matrix of i,j points E, use scale() to center, then var/covar matrix F=E'E--eigenanalyze for second moment of area

# spatial eigenfunction analysis
# to troubleshoot: have to use function from "ff" package to save distance matrices to disk--they're going to be huge (>150GB)
gen.eignf<-function(x,cortical=FALSE) {
  require(ape)
  require(vegan)
  require(ff)
  if(cortical) {
    posdist<-dist(x,method="manhattan")
    posdist[which(posdist>1)]<-4
    sp.eigenf<-pcoa(posdist)    
  } else {
    posdist<-dist(x,method="euclidean")
    postree<-spantree(posdist)
    posmax<-max(postree$dist)
    posdist[which(posdist>posmax)]<-posmax*4
    sp.eigenf<-pcoa(posdist)  
  }
  result<-sp.eigenf$vectors
  return(result)
}


## Create a formula for a model with a large number of variables:
xnam <- paste0("x", 1:25)
(fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"))))

# test section Gavia immer ulna 50.16um thick, 546 nm, birefringence 0.005
Gavia<-read.Rotopol.spec(50.16,546,0.005, mask=TRUE)
