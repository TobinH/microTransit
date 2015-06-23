# read Rotopol (Werner Kaminsky) bitmaps into a qPLM object in R 
# with specified parameters.
# arguments--

# sample.name: string for file name of exported images

# bitmap.blob: string to direct Sys.glob() to I, |sind|, and phi bitmaps.

# mask.glob: string to direct Sys.glob() to mask bitmap.

# nw.thickness - ne.thickness: slide thickness in microns at top left (NW),
# bottom right (SE), bottom left (SW), and top right (NE) of masked element.

# wavelength: center of analysis wavelength in nm. Default is 532 nm.

# birefringence: estimated birefringence constant of ROI. Default is 0.005
#    (empirically determined for extant bone tissue).

# pixel: pixel size in microns--default is for E145 Leica A6 scope at 1.25x

# up: if the image is skewed, up is whatever direction you want to assign 
#   to the top, as a compass heading in degrees
#   (1-359--extremes are 90 for the right side, 180 for bottom, 270 for left)

# mask: Loads a black-and-white bitmap that defines a region of interest. 
#   Black pixels (0) will be dropped, 
#   to allow you to choose the mask file.

# pics.out: if TRUE, returns single-channel greyscale .tiffs 
#   and a PolarLUV composite .tiff

# theta.scale: if TRUE, scales "_overview.tif" image color scale to
#   peak at maximum theta value (default is FALSE).

# theta.max: if theta.scale is FALSE, specified upper end (degrees) 
#   for "_overview.tif" color scale (default is 90 degrees).

# invert.theta.scale: if TRUE, bright colors pointing toward viewer, 
#   brightness falls off as orientations go to parallel with the slide. 
#   Useful for sections with predominantly out-of-plane fibers 
#   (default is FALSE).

Rotopol.qPLM.script.CMCs<-function(sample.name=NULL,
                       bitmap.glob,
                       mask.glob=NULL,
                       nw.thickness, 
                       se.thickness,
                       sw.thickness,
                       ne.thickness,
                       wavelength=532, 
                       birefringence=0.005, 
                       pixel=7.5832259, 
                       up=0, 
                       mask=TRUE, 
                       pics.out=TRUE,
                       theta.scale=FALSE,
                       theta.max=90,
                       invert.theta.scale=FALSE) {
  require(bmp)
  require(colorspace)
  require(EBImage)
  
  Rotopol.raw<-vector("list", 7)
  # create an object to hold Rotopol bitmap values 
  # and other measured/estimated parameters
  
  bmpfiles<-Sys.glob(bitmap.glob)
  # Rotopol bitmap file selection 
  
  if (mask) {
    maskfile<-Sys.glob(mask.glob)
  }
  # mask bitmap file selection
  
  bmpfiles<-bmpfiles[order(bmpfiles)]
  # slots Rotopol a0-a2 bitmaps into numerical order
  
  pb<-winProgressBar(title="Rotopol.PLM", min = 0,
                       max = 19, width = 300)
  # (Windows) create a progress bar--
  # mostly an instant-feedback debugging tool, 
  # some operations can be time-consuming
  
  setWinProgressBar(pb, 1, title="Rotopol.qPLM: reading I")
  Rotopol.raw[[1]]<-read.bmp(bmpfiles[1])
  # read a0 bitmap (transmittance, I)
  
  setWinProgressBar(pb, 2, title="Rotopol.qPLM: reading |sin d|")
  Rotopol.raw[[2]]<-read.bmp(bmpfiles[2])
  # read a1 bitmap (retardance, |sin d|)
  
  setWinProgressBar(pb, 3, title="Rotopol.qPLM: reading phi")
  Rotopol.raw[[3]]<-read.bmp(bmpfiles[3])
  # read a2 bitmap (azimuth, phi)
  
  setWinProgressBar(pb, 4, title="Rotopol.qPLM: applying mask")
      
  ifelse(mask, 
         Rotopol.raw[[7]]<-read.bmp(maskfile)/255, 
         Rotopol.raw[[7]]<-array(1,dim=dim(Rotopol.raw[[1]])))
  # read mask or create blank mask object
    
  dindx<-which(Rotopol.raw[[7]][,,2]>1/255, arr.ind=TRUE)
  # index of masked pixels
  
  ifelse(as.logical(min(dindx[,1])<10),
         umin<-1,
         umin<-min(dindx[,1])-10)
  ifelse(as.logical(max(dindx[,1])>dim(Rotopol.raw[[7]])[1]-10),
         umax<-dim(Rotopol.raw[[7]])[1],
         umax<-max(dindx[,1])+10)
  ifelse(as.logical(min(dindx[,2])<10),
         vmin<-1,
         vmin<-min(dindx[,2])-10)
  ifelse(as.logical(max(dindx[,2])>dim(Rotopol.raw[[7]])[2]-10),
         vmax<-dim(Rotopol.raw[[7]])[2],
         vmax<-max(dindx[,2])+10)
  
  Rotopol.raw[[7]]<-as.matrix(Rotopol.raw[[7]][umin:umax, vmin:vmax,2])
  Rotopol.raw[[1]]<-as.matrix(Rotopol.raw[[1]][umin:umax, vmin:vmax,2]
                              *Rotopol.raw[[7]])
  Rotopol.raw[[2]]<-as.matrix(Rotopol.raw[[2]][umin:umax, vmin:vmax,2]
                              *Rotopol.raw[[7]])
  Rotopol.raw[[3]]<-as.matrix(Rotopol.raw[[3]][umin:umax, vmin:vmax,2]
                              *Rotopol.raw[[7]])
  # trim image arrays to non-zero pixels with a 10-pixel blank border
    
  Rotopol.distilled<-array(data=NA, 
                           c(dim(Rotopol.raw[[1]])[2], 
                             dim(Rotopol.raw[[1]])[1], 3))
  # create output object array
  
  setWinProgressBar(pb, 5, title="Rotopol.qPLM: applying parameters")
  
  thickness<-as.numeric(c(nw.thickness,
                          se.thickness,
                          sw.thickness,
                          ne.thickness))
  Rotopol.raw[[5]]<-as.numeric(wavelength)
  Rotopol.raw[[6]]<-as.numeric(birefringence)
  # parameters from arguments
  
   wedge<-matrix(c(1, 1,
                  umax-umin+1, vmax-vmin+1,
                  umax-umin+1, 1,
                  1, vmax-vmin+1), 
                byrow=TRUE, nrow=4, 
                dimnames=list(c("NW","SE", "SW", "NE"), c("u", "v")))
  # NW, SE, SW, NE positions in (u, v) coordinates
  
  wedge<-as.data.frame(cbind(thickness, wedge))
  # data frame for determining wedging
  
  wedge.f<-lm(thickness ~ u+v, data=wedge)
  # linear fit for wedging
  
  u.v.pos<-expand.grid(1:(umax-umin+1), 1:(vmax-vmin+1))
  colnames(u.v.pos)<-c("u", "v")
  # pixel positions for by-pixel thickness estimate
  
  Rotopol.raw[[4]]<-matrix(predict(wedge.f, u.v.pos), nrow=umax-umin+1)
  # by-pixel thickness estimates
  
  setWinProgressBar(pb, 6, title="Rotopol.qPLM: writing array: I")
  Rotopol.distilled[,,1]<-t(Rotopol.raw[[1]]/255)
  # transmittance pixels scaled to 0-1 range
  
  setWinProgressBar(pb, 7, title="Rotopol.qPLM: writing array: Theta")
  Rotopol.distilled[,,2]<-t(asin
                             (sqrt
                                  (
                                    (asin(Rotopol.raw[[2]]/255))
                                    *Rotopol.raw[[5]]
                                    /(2*pi*Rotopol.raw[[4]]
                                      *1000
                                      *Rotopol.raw[[6]])
                                    )))
  # retardance pixels transformed to elevation angle (theta) 
  # mapped linearly to 0-1 range
  
  setWinProgressBar(pb, 8, title="Rotopol.qPLM: writing array: Phi")
  Rotopol.distilled[,,3]<-t(Rotopol.raw[[3]]/255)
  # azimuth pixels (0 degrees to 179 degrees) scaled to 0-1 range

  
  attributes(Rotopol.distilled)<-list(dim=dim(Rotopol.distilled),
                                      dimnames=list(paste("u",seq(length.out=length(Rotopol.distilled[,1,1])), sep=""), 
                                                    paste("v", seq(length.out=length(Rotopol.distilled[1,,1])), sep=""),
                                                    c("Trans","Theta","Phi")),
                                      thickness_um=Rotopol.raw[[4]],
                                      wavelength_nm=wavelength,
                                      birefringence=birefringence,
                                      pixel.size_um=pixel,
                                      ccw.skew_deg=up)
  # sets attributes for qPLM object
 
  rm(Rotopol.raw)
  gc()
  
  if (pics.out) {
    setWinProgressBar(pb, 10, title="Rotopol.qPLM: building colorspace")
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
    
    setWinProgressBar(pb, 11, title="Rotopol.qPLM: LUV encoding")
    PLUV.LUT<-polarLUV((pixLUT[,1]*0.75+22.5)*90/theta.max, 
                       pixLUT[,1]*57.65/theta.max, 
                       pixLUT[,2]*360/180)
    # polar LUV colorspace encoding of each integer combination
    
    setWinProgressBar(pb, 12, title="Rotopol.qPLM: RGB translation")
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
    
    setWinProgressBar(pb, 13, title="Rotopol.qPLM: matching orientation to color")
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
    
    coded<-array(data=RGBmat[LUindex,], dim=dim(Rotopol.distilled))
    rm(RGBmat)
    rm(LUindex)
    gc()
    # cast looked-up RGB values into RGB image array
        
    setWinProgressBar(pb, 14, title="Rotopol.qPLM: writing image: _trans")
    trans<-Image(Rotopol.distilled[,,1], colormode="Grayscale")
    writeImage(trans, file=paste(sample.name,"_trans.tif", sep=""), bits.per.sample=8L, type="tiff")
    # Transmission (I) grayscale image out
    
    setWinProgressBar(pb, 15, title="Rotopol.qPLM: writing image: _theta")
    theta<-Image(Rotopol.distilled[,,2], colormode="Grayscale")
    writeImage(theta, file=paste(sample.name, "_theta.tif", sep=""), bits.per.sample=8L, type="tiff")
    # colatitude (theta) grayscale image out
    
    setWinProgressBar(pb, 16, title="Rotopol.qPLM: writing image: _phi")
    phi<-Image(Rotopol.distilled[,,3], colormode="Grayscale")
    writeImage(phi, file=paste(sample.name, "_phi.tif",sep=""), bits.per.sample=8L, type="tiff")
    rm(trans, theta, phi)
    gc()
    # azimuth (phi) grayscale image out
    
    setWinProgressBar(pb, 17, title="Rotopol.qPLM: writing image: _overview")
    output<-Image(coded,colormode="Color")
    writeImage(output,file=paste(sample.name,"_overview.tif", sep=""), bits.per.sample=8L, type="tiff")
    # combined colatitude and azimuth color image out
    
    setWinProgressBar(pb, 18, title="Rotopol.qPLM: showing overview")
    display(output)
    # pull up overview image in browser window
    
  }
  setWinProgressBar(pb, 19, title="Rotopol.qPLM: done!")
  close(pb)
  
  invisible(Rotopol.distilled)
  return(Rotopol.distilled)
  # sends calibrated qPLM data to specified object
  
  gc()
  # may be overkill, but see above
}

# test<-Rotopol.qPLM("test", 15, 15, 15, 15, 532, 0.005)