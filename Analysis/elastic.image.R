# Elastic model plots through EBImage (probably much faster)

fig.elastic<-function(elastic,
                      sample.name,
                      pixel=7.5832259){
  require(EBImage)
  elastic[,1]<-elastic[,1]+abs(min(elastic[,1]))
  elastic[,2]<-elastic[,2]+abs(min(elastic[,2]))
  elastic[,1:2]<-elastic[,1:2]/pixel
  elastic[,1:2]<-elastic[,1:2]*1000
  #return to u,v pixel scaling
  
  scale<-max(abs(elastic[,3:5]))
  elastic[,3:5]<-ceiling(((elastic[,3:5]/max(abs(elastic[,3:5])))+1)*127.5)
  # scale strains to grayscale output
  
  elastic.images<-array(data=0, dim=c(max(elastic[,2], max(elastic[,1]),3)))
  print(length(elastic.images)[1])
  print(length(elastic.images)[2])
  # setup results array
  
  for (i in 1:length(elastic.images)[1]){
    for (j in 1:length(elastic.images)[2]){
      elastic.images[i,j,]<-elastic[which(i %in% elastic[,2] & j %in% elastic[,1],3:5)]
    }
  }
  
  iso<-Image(elastic.images[,,1])
  lo<-Image(elastic.images[,,2])
  diff<-Image(elastic.images[,,3])
  writeImage(iso, file=paste(sample.name, "_elastic_iso.tif",sep=""), bits.per.sample=8L, type="tiff")
  writeImage(lo, file=paste(sample.name, "_elastic_lo.tif",sep=""), bits.per.sample=8L, type="tiff")
  writeImage(diff, file=paste(sample.name, "_elastic_diff.tif",sep=""), bits.per.sample=8L, type="tiff")
  display(iso)
  display(lo)
  display(diff)
}