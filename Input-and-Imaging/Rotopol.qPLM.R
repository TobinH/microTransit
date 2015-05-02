# read Rotopol (Werner Kaminsky) bitmaps into a qPLM object in R 
# with specified parameters.
# arguments--
# thickness: slide thickness in microns
# wavelength: center of analysis wavelength in nm
# birefringence: estimated birefringence constant of ROI
# pixel: pixel size in microns--default is for E145 Leica A6 scope at 1.25x
# up: if the image is skewed, up is whatever direction you want to assign 
#   to the top, as a compass heading in degrees
#   (1-359--extremes are 90 for the right side, 180 for bottom, 270 for left)
# mask: Loads a black-and-white bitmap that defines a region of interest. 
#   Black pixels (0) will be dropped, 
#   to allow you to choose the mask file.
# pics.out: if TRUE, returns single-channel greyscale .tiffs 
#   and a PolarLUV composite .tiff
# mono.bkgrnd: Gray level for single channel backgrounds (0-100). 
#   Default is white (0).
# comp.bkgrnd: Gray level for composite background (0-100, 100 is black). Default 
#   is 50% gray (50)
# theta.scale: if TRUE, scales "_overview.tif" image color scale to  peak at maximum theta value
#   (default is FALSE)
# theta.max: if theta.scale is FALSE, specified upper end (degrees) for "_overview.tif" color scale
#   (default is 90 degrees).
# invert.theta.scale: if TRUE, bright colors pointing toward viewer, brightness falls off as orientations go
#   to parallel with the slide. Useful for sections with predominantly out-of-plane fibers (default is FALSE)

Rotopol.qPLM<-function(sample.name, 
                       thickness, 
                       wavelength, 
                       birefringence, 
                       pixel=7.5832259, 
                       up=0, 
                       mask=FALSE, 
                       pics.out=TRUE, 
                       comp.bkgrnd=50,
                       theta.scale=FALSE,
                       theta.max=90,
                       invert.theta.scale=FALSE) {
  require(bmp)
  require(colorspace)
  require(EBImage)
  
  Rotopol.raw<-vector("list", 7)
  # create an object to hold Rotopol bitmap values and other measured/estimated parameters
  
  bmpfiles<-choose.files(default="", caption="Select a0, a1, and a2 bitmaps",
                         multi=TRUE)
  # Rotopol bitmap file selection dialog
  
  if (mask) {
    maskfile<-choose.files(default="", caption="Select mask bitmap")
  }
  # mask bitmap file selection dialog
  
  bmpfiles<-bmpfiles[order(bmpfiles)]
  # slots Rotopol a0-a2 bitmaps into numerical order
  
  pb <- winProgressBar(title="Rotopol.PLM", min = 0,
                       max = 19, width = 300)
  # (Windows) create a progress bar--mostly an instant-feedback debugging tool, some operations can be time-consuming
  
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
  # parameters from arguments
  
  ifelse(mask, Rotopol.raw[[7]]<-read.bmp(maskfile)/255, Rotopol.raw[[7]]<-array(1,dim=dim(Rotopol.raw[[1]])))
  # read mask or create blank mask object
  
  Rotopol.distilled<-array(data=NA, c(dim(Rotopol.raw[[1]])[2],dim(Rotopol.raw[[1]])[1],3))
  # create ouput object array
  
  setWinProgressBar(pb, 5, title="Rotopol.qPLM: writing array: I")
  Rotopol.distilled[,,1]<-t(as.matrix(Rotopol.raw[[1]][,,2]/255))
  # transmittance pixels scaled to 0-1 range
  
  setWinProgressBar(pb, 6, title="Rotopol.qPLM: writing array: Theta")
  Rotopol.distilled[,,2]<-t(as.matrix(asin(sqrt((asin(Rotopol.raw[[2]][,,2]/255))
                                                *Rotopol.raw[[5]]/(2*pi*Rotopol.raw[[4]]*1000*Rotopol.raw[[6]])))))
  # retardance pixels transformed to elevation angle (theta) mapped linearly to 0-1 range
  
  setWinProgressBar(pb, 7, title="Rotopol.qPLM: writing array: Phi")
  Rotopol.distilled[,,3]<-t(as.matrix(Rotopol.raw[[3]][,,2]/255))
  # azimuth pixels (0 degrees to 179 degrees) scaled to 0-1 range
  
  setWinProgressBar(pb, 8, title="Rotopol.qPLM: applying mask")
  Rotopol.distilled[,,1]<-t(Rotopol.raw[[7]][,,2])*Rotopol.distilled[,,1]
  Rotopol.distilled[,,2]<-t(Rotopol.raw[[7]][,,2])*Rotopol.distilled[,,2]
  Rotopol.distilled[,,3]<-t(Rotopol.raw[[7]][,,2])*Rotopol.distilled[,,3]
  # mask cutouts--only white mask pixels remain
  
  rm(Rotopol.raw)
  gc()
  # manipulating these arrays can tax some 32-bit systems... better safe and slow than
  # lost in the woods with a programmer that doesn't have good error handling skills (me)
  
  attributes(Rotopol.distilled)<-list(dim=dim(Rotopol.distilled),
                                      dimnames=list(paste("v",seq(length.out=length(Rotopol.distilled[,1,1])), sep=""), 
                                                    paste("w", seq(length.out=length(Rotopol.distilled[1,,1])), sep=""),
                                                    c("Trans","Theta","Phi")),
                                      thickness_um=thickness,
                                      wavelength_nm=wavelength,
                                      birefringence=birefringence,
                                      pixel.size_um=pixel,
                                      ccw.skew_deg=up)
  # sets attributes for qPLM object
  
  if (pics.out) {
    setWinProgressBar(pb, 9, title="Rotopol.qPLM: building colorspace")
    if(theta.scale){
      theta.max<-ceiling(max(Rotopol.distilled[,,2]*90))
    }
    theta.max<-ceiling(theta.max)
    if (invert.theta.scale) {
      thetaseq<-seq(from=90, to=90-theta.max)
    }
    else{
      thetaseq<-seq(from=0, to=theta.max) 
    }
    phiseq<-seq(from=0, to=180)
    pixLUT<-expand.grid(thetaseq, phiseq)
    rm(thetaseq, phiseq)
    # integer combinations of theta and phi for the look-up table
    
    setWinProgressBar(pb, 10, title="Rotopol.qPLM: LUV encoding")
    PLUV.LUT<-polarLUV(pixLUT[,1]*73.2/theta.max, pixLUT[,1]*57.65/theta.max, pixLUT[,2]*360/180)
    # polar LUV colorspace encoding of each integer combination
    
    setWinProgressBar(pb, 11, title="Rotopol.qPLM: RGB translation")
    RGB.LUT<-as(PLUV.LUT,"RGB")
    rm(PLUV.LUT)
    # RGB value conversion of polar LUV look-up table
    
    pixmat<-as.matrix(pixLUT)
    rm(pixLUT)
    # n by 2 integer combination table
    
    RGBmat<-RGB.LUT@coords
    rm(RGB.LUT)
    gc()
    # n by 3 RGB look-up table
    
    setWinProgressBar(pb, 12, title="Rotopol.qPLM: matching orientation to color")
    if (invert.theta.scale){
      theta.values<-as.vector(as.integer((1-Rotopol.distilled[,,2])*90))
    }
    else{
      theta.values<-as.vector(as.integer(Rotopol.distilled[,,2]*90))
    }
      
    immat<-cbind(theta.values, as.vector((as.integer(Rotopol.distilled[,,3]*180))))
    # image pixel values m by 2 table
    
    LUindex<-match(data.frame(t(immat)), data.frame(t(pixmat)))
    rm(immat, pixmat)
    # look-up indices for image pixels
    
    masked<-array(data=RGBmat[LUindex,], dim=dim(Rotopol.distilled))
    rm(RGBmat)
    rm(LUindex)
    gc()
    # cast looked-up RGB values into RGB image array
    
    setWinProgressBar(pb, 13, title="Rotopol.qPLM: trimming images")
    dindx<-which(masked>1/256, arr.ind=TRUE)
    dindx<-unique(dindx[,1:2])
    # index of non-zero pixel values for theta and phi
    
    ifelse(as.logical(min(dindx[,1])<10),
           wmin<-0,
           wmin<-min(dindx[,1]))
    ifelse(as.logical(max(dindx[,1])>dim(Rotopol.distilled)[1]-10),
           wmax<-dim(Rotopol.distilled)[1],
           wmax<-max(dindx[,1]))
    ifelse(as.logical(min(dindx[,2])<10),
           vmin<-0,
           vmin<-min(dindx[,2]))
    ifelse(as.logical(max(dindx[,2])>dim(Rotopol.distilled)[2]-10),
           vmax<-dim(Rotopol.distilled)[2],
           vmax<-max(dindx[,2]))
    masked<-masked[wmin:wmax, vmin:vmax,]
    # trim image array to non-zero pixels with a 10-pixel blank border
    
    setWinProgressBar(pb, 14, title="Rotopol.qPLM: writing image: _trans")
    trans<-Image(Rotopol.distilled[wmin:wmax,vmin:vmax,1], colormode="Grayscale")
    writeImage(trans, file=paste(sample.name,"_trans.tif", sep=""), bits.per.sample=8L, type="tiff")
    # Transmission (I) grayscale image out
    
    setWinProgressBar(pb, 15, title="Rotopol.qPLM: writing image: _theta")
    theta<-Image(Rotopol.distilled[wmin:wmax, vmin:vmax,2], colormode="Grayscale")
    writeImage(theta, file=paste(sample.name, "_theta.tif", sep=""), bits.per.sample=8L, type="tiff")
    # colatitude (theta) grayscale image out
    
    setWinProgressBar(pb, 16, title="Rotopol.qPLM: writing image: _phi")
    phi<-Image(Rotopol.distilled[wmin:wmax,vmin:vmax,3], colormode="Grayscale")
    writeImage(phi, file=paste(sample.name, "_phi.tif",sep=""), bits.per.sample=8L, type="tiff")
    rm(trans, theta, phi)
    gc()
    # azimuth (phi) grayscale image out
    
    setWinProgressBar(pb, 17, title="Rotopol.qPLM: writing image: _overview")
    out<-array(data=(100-comp.bkgrnd)/100, dim=dim(masked))
    # create a blank grey background (adjustable in args)
    
    for (i in 1:nrow(masked)){
      dindx<-which(masked[i,,]>1/256, arr.ind=TRUE)
      dindx<-unique(dindx[,1])
      out[i,dindx,]<-masked[i,dindx,]
    }  
    rm(masked)
    gc()
    # overlaying RGB theta and phi combined image on null pixel background
    
    output<-Image(out,colormode="Color")
    writeImage(output,file=paste(sample.name,"_overview.tif", sep=""), bits.per.sample=8L, type="tiff")
    # combined colatitude and azimuth color image out
    
    setWinProgressBar(pb, 18, title="Rotopol.qPLM: showing overview")
    display(output)
    # pull up overview image in browser window
    
    print("calibrated images written to working directory")
  }
  setWinProgressBar(pb, 19, title="Rotopol.qPLM: done!")
  
  invisible(Rotopol.distilled)
  return(Rotopol.distilled)
  # sends calibrated qPLM data to specified object
  
  close(pb)
  gc()
  # may be overkill, but see above
}

# test<-Rotopol.qPLM("test", 15, 532, 0.005)