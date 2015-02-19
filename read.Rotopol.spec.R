# read rotopol into R with specified parameters

read.Rotopol.spec<-function(thickness,wavelength,birefringence,mask=FALSE) {
  require(bmp)
  Rotopol.raw<-vector("list", 7)
  bmpfiles<-choose.files(default="",caption="Select a0, a1, and a2 bitmaps",multi=TRUE)
  if (mask) {
    maskfile<-choose.files(default="",caption="Select mask bitmap")
  }
  bmpfiles<-bmpfiles[order(bmpfiles)]
  Rotopol.raw[[1]]<-read.bmp(bmpfiles[1])
  Rotopol.raw[[2]]<-read.bmp(bmpfiles[2])
  Rotopol.raw[[3]]<-read.bmp(bmpfiles[3])
  Rotopol.raw[[4]]<-as.numeric(thickness)
  Rotopol.raw[[5]]<-as.numeric(wavelength)
  Rotopol.raw[[6]]<-as.numeric(birefringence)
  ifelse(mask, Rotopol.raw[[7]]<-read.bmp(maskfile)/256, Rotopol.raw[[7]]<-array(1,dim=dim(Rotopol.raw[[1]])))
  Rotopol.distilled<-array(data=NA, c(dim(Rotopol.raw[[1]])[2],dim(Rotopol.raw[[1]])[1],3))
  Rotopol.distilled[,,1]<-t(as.matrix(Rotopol.raw[[1]][,,2]/256))
  Rotopol.distilled[,,2]<-t(as.matrix(sqrt(asin(Rotopol.raw[[2]][,,2]/256)*(Rotopol.raw[[5]]/(2*pi*Rotopol.raw[[4]]*1000*Rotopol.raw[[6]])))))
  Rotopol.distilled[,,3]<-t(as.matrix(Rotopol.raw[[3]][,,2]/256))
  Rotopol.distilled[,,1]<-t(Rotopol.raw[[7]][,,2])*Rotopol.distilled[,,1]
  Rotopol.distilled[,,2]<-t(Rotopol.raw[[7]][,,2])*Rotopol.distilled[,,2]
  Rotopol.distilled[,,3]<-t(Rotopol.raw[[7]][,,2])*Rotopol.distilled[,,3]
  return(Rotopol.distilled)
}

# test<-read.Rotopol.spec(15,546,0.005)
