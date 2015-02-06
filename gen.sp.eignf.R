gen.sp.eignf<-function(x,cortical=FALSE) {
  require(ape)
  require(vegan)
  if(cortical) {
    posdist<-dist(x,method="manhattan")
    print(paste("Distance matrix completed at",Sys.time()))
    posdist[which(posdist>1)]<-4
    sp.eigenf<-pcoa(posdist)
    print(paste("PCoA completed at",Sys.time()))
  } else {
    posdist<-dist(x,method="euclidean")
    print(paste("Distance matrix completed at",Sys.time()))
    postree<-spantree(posdist)
    print(paste("Minimum spanning tree completed at",Sys.time()))
    posmax<-max(postree$dist)
    posdist[which(posdist>posmax)]<-posmax*4
    sp.eigenf<-pcoa(posdist)
    print(paste("PCoA completed at",ys.time()))
  }
