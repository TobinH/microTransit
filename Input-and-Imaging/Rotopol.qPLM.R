# read Rotopol (Werner Kaminsky) bitmaps into a qPLM object in R with specified parameters
# arguments--
# thickness: slide thickness in microns
# wavelength: center of analysis wavelength in nm
# birefringence: estimated birefringence constant of ROI
# pixel: pixel size in microns--default is for E145 Leica A6 scope at 1.25x
# up: if the image is skewed, up is whatever direction you want to assign to the top, as a compass heading in degrees
#   (1-359--extremes are 90 for the right side, 180 for bottom, 270 for left)
# mask: Loads a black-and-white bitmap that defines a region of interest. Black pixels (0) will be dropped, 
#   white pixels (255) will be kept. If mask=TRUE, a window will open to allow you to choose the mask file.
# pics.out: if TRUE, returns single-channel greyscale .tiffs and a PolarLUV composite .tiff
# mono.bkgrnd: Gray level for single channel backgrounds (0-100). Default is white (0).
# comp.bkgrnd: Gray level for composite background (0-100). Default is 50% gray (50)

Rotopol.qPLM<-function(sample.name,thickness, wavelength, birefringence, pixel=7.5832259, up=0, mask=FALSE, pics.out=TRUE, comp.bkgrnd=50) {
  require(bmp)
  require(colorspace)
  require(EBImage)
  Rotopol.raw<-vector("list", 7)
  bmpfiles<-choose.files(default="",caption="Select a0, a1, and a2 bitmaps",multi=TRUE)
  if (mask) {
    maskfile<-choose.files(default="",caption="Select mask bitmap")
  }
  bmpfiles<-bmpfiles[order(bmpfiles)]
  pb <- winProgressBar(title="Rotopol.PLM", min = 0,
                       max = 19, width = 300)
  setWinProgressBar(pb, 1, title="Rotopol.qPLM: reading I")
  Rotopol.raw[[1]]<-read.bmp(bmpfiles[1])
  # read a0 bitmap (transmittance, I)
  setWinProgressBar(pb, 2, title="Rotopol.qPLM: reading |sin d|")
  Rotopol.raw[[2]]<-read.bmp(bmpfiles[2])
  # read a1 bitmap (retardance, |sin d|)
  setWinProgressBar(pb, 3, title="Rotopol.qPLM: reading phi")
  Rotopol.raw[[3]]<-read.bmp(bmpfiles[3])
  # read a2 bitmap (azimuth, phi)
  setWinProgressBar(pb, 4, title="Rotopol.qPLM: parameters & mask")
  Rotopol.raw[[4]]<-as.numeric(thickness)
  Rotopol.raw[[5]]<-as.numeric(wavelength)
  Rotopol.raw[[6]]<-as.numeric(birefringence)
  ifelse(mask, Rotopol.raw[[7]]<-read.bmp(maskfile)/256, Rotopol.raw[[7]]<-array(1,dim=dim(Rotopol.raw[[1]])))
  Rotopol.distilled<-array(data=NA, c(dim(Rotopol.raw[[1]])[2],dim(Rotopol.raw[[1]])[1],3))
  setWinProgressBar(pb, 5, title="Rotopol.qPLM: writing array: I")
  Rotopol.distilled[,,1]<-t(as.matrix(Rotopol.raw[[1]][,,2]/255))
  # transmittance pixels scaled to 0-1 range
  setWinProgressBar(pb, 6, title="Rotopol.qPLM: writing array: Theta")
  Rotopol.distilled[,,2]<-t(as.matrix(sqrt(asin(Rotopol.raw[[2]][,,2]/255)*(Rotopol.raw[[5]]/(2*pi*Rotopol.raw[[4]]*1000*Rotopol.raw[[6]])))))
  # retardance pixels transformed to elevation angle (theta)-theta mapped linearly to 0-1 range
  setWinProgressBar(pb, 7, title="Rotopol.writing array: Phi")
  Rotopol.distilled[,,3]<-t(as.matrix(Rotopol.raw[[3]][,,2]/255))
  # azimuth pixels (0 degrees to 179 degrees) scaled to 0-1 range
  setWinProgressBar(pb, 8, title="Rotopol.qPLM: applying mask")
  Rotopol.distilled[,,1]<-t(Rotopol.raw[[7]][,,2])*Rotopol.distilled[,,1]
  Rotopol.distilled[,,2]<-t(Rotopol.raw[[7]][,,2])*Rotopol.distilled[,,2]
  Rotopol.distilled[,,3]<-t(Rotopol.raw[[7]][,,2])*Rotopol.distilled[,,3]
  # mask cutouts--only white mask pixels remain
  rm(Rotopol.raw)
  attributes(Rotopol.distilled)<-list(dim=dim(Rotopol.distilled),dimnames=list(paste("v",seq(length.out=length(Rotopol.distilled[,1,1])),sep=""),paste("w",seq(length.out=length(Rotopol.distilled[1,,1])),sep=""),c("Trans","Theta","Phi")),thickness_um=thickness,wavelength_nm=wavelength,birefringence=birefringence,pixel.size_um=pixel,ccw.skew_deg=up)
  # sets attributes for qPLM object
  if (pics.out) {
    setWinProgressBar(pb, 9, title="Rotopol.qPLM: building colorspace")
    thetaseq<-seq(from=0, to=90)
    phiseq<-seq(from=0, to=180)
    pixLUT<-expand.grid(thetaseq,phiseq)
    rm(thetaseq,phiseq)
    setWinProgressBar(pb, 10, title="Rotopol.qPLM: LUV encoding")
    PLUV.LUT<-polarLUV(pixLUT[,1]*73.2/90,pixLUT[,1]*57.65/90,pixLUT[,2]*360/180)
    setWinProgressBar(pb, 11, title="Rotopol.qPLM: RGB translation")
    RGB.LUT<-as(PLUV.LUT,"RGB")
    rm(PLUV.LUT)
    pixmat<-as.matrix(pixLUT)
    rm(pixLUT)
    RGBmat<-RGB.LUT@coords
    rm(RGB.LUT)
    gc()
    setWinProgressBar(pb, 12, title="Rotopol.qPLM: matching orientation to color")
    immat<-cbind(as.vector(as.integer(Rotopol.distilled[,,2]*90)),as.vector((as.integer(Rotopol.distilled[,,3]*180))))
    LUindex<-match(data.frame(t(immat)),data.frame(t(pixmat)))
    rm(immat,pixmat)
    masked<-array(data=RGBmat[LUindex,],dim=dim(Rotopol.distilled))
    rm(RGBmat)
    rm(LUindex)
    gc()
    setWinProgressBar(pb, 13, title="Rotopol.qPLM: trimming images")
    dindx<-which(masked>1/256,arr.ind=TRUE)
    dindx<-unique(dindx[,1:2])
    ifelse(as.logical(min(dindx[,1])<10),wmin<-0,wmin<-min(dindx[,1]))
    ifelse(as.logical(max(dindx[,1])>dim(Rotopol.distilled)[1]-10),wmax<-dim(Rotopol.distilled)[1],wmax<-max(dindx[,1]))
    ifelse(as.logical(min(dindx[,2])<10),vmin<-0,vmin<-min(dindx[,2]))
    ifelse(as.logical(max(dindx[,2])>dim(Rotopol.distilled)[2]-10),vmax<-dim(Rotopol.distilled)[2],vmax<-max(dindx[,2]))
    masked<-masked[wmin:wmax,vmin:vmax,]
    setWinProgressBar(pb, 14, title="Rotopol.qPLM: writing image: _trans")
    trans<-Image(Rotopol.distilled[wmin:wmax,vmin:vmax,1],colormode="Grayscale")
    writeImage(trans, file=paste(sample.name,"_trans.tif",sep=""), type="tiff")
    setWinProgressBar(pb, 15, title="Rotopol.qPLM: writing image: _theta")
    theta<-Image(Rotopol.distilled[wmin:wmax,vmin:vmax,2],colormode="Grayscale")
    writeImage(theta, file=paste(sample.name,"_theta.tif",sep=""), type="tiff")
    setWinProgressBar(pb, 16, title="Rotopol.qPLM: writing image: _phi")
    phi<-Image(Rotopol.distilled[wmin:wmax,vmin:vmax,3],colormode="Grayscale")
    writeImage(phi, file=paste(sample.name,"_phi.tif",sep=""), type="tiff")
    rm(trans,theta,phi)
    gc()
    setWinProgressBar(pb, 17, title="Rotopol.qPLM: writing image: _overview")
    out<-array(data=comp.bkgrnd/100,dim=dim(masked))
    for (i in 1:nrow(masked)){
      dindx<-which(masked[i,,]>1/256,arr.ind=TRUE)
      dindx<-unique(dindx[,1])
      out[i,dindx,]<-masked[i,dindx,]
    }
    rm(masked)
    gc()
    output<-Image(out,colormode="Color")
    writeImage(output,file=paste(sample.name,"_overview.tif",sep=""), bits.per.sample=8L, type="tiff")
    setWinProgressBar(pb, 18, title="Rotopol.qPLM: showing overview")
    display(output)
    setWinProgressBar(pb, 19, title="Rotopol.qPLM: done!")
    print("calibrated images written to working directory")
    close(pb)
  }
  invisible(Rotopol.distilled)
  return(Rotopol.distilled)
}

# test<-read.Rotopol.spec(15,546,0.005)
