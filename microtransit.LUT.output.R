# polar coord palette generation to create a LUT for microtransit output
# will hopefully make composite images faster

microtransit.LUT.output<-function(x,sample.name){
  require(colorspace)
  require(EBImage)
  #trans<-Image(x[,,1],colormode="Grayscale")
  #theta<-Image(x[,,2],colormode="Grayscale")
  #phi<-Image(x[,,3],colormode="Grayscale")
  #writeImage(trans, file=paste(sample.name,"_trans.tif",sep=""), type="tiff")
  #writeImage(theta, file=paste(sample.name,"_theta.tif",sep=""), type="tiff")
  #writeImage(phi, file=paste(sample.name,"_phi.tif",sep=""), type="tiff")
  thetaseq<-seq(from=0, to=90)
  phiseq<-seq(from=0, to=180)
  pixLUT<-expand.grid(thetaseq,phiseq)
  rm(thetaseq)
  rm(phiseq)
  PLUV.LUT<-polarLUV(pixLUT[,1]*73.2/90,pixLUT[,1]*57.65/90,pixLUT[,2]*360/180)
  RGB.LUT<-as(PLUV.LUT,"RGB")
  rm(PLUV.LUT)
  pixmat<-as.matrix(pixLUT)
  rm(pixLUT)
  RGBmat<-RGB.LUT@coords
  rm(RGB.LUT)
  gc()
  immat<-cbind(as.vector(as.integer(x[,,2]*90)),as.vector((as.integer(x[,,3]*180))))
  LUindex<-match(data.frame(t(immat)),data.frame(t(pixmat)))
  output<-array(data=RGBmat[LUindex,],dim=dim(x))
  bkgrnd<-array(data=0.5,dim=dim(x))
  bkgrnd[which(output>5/256)]<-output[which(output>5/256)]
  output<-bkgrnd
  rm(bkgrnd)
  rm(immat)
  rm(LUindex)
  gc()
  output<-Image(output,colormode="Color")
  gc()
  #display(outimg)
  writeImage(output,file=paste(sample.name,"_overview.tif",sep=""), bits.per.sample=8L, type="tiff")
  #print("calibrated images written to working directory")
  return(output)
}
