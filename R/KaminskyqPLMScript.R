#' @title Build qPLM Array Object From Specified Rotopol(tm) Format Images
#'
#' @description
#' \code{KaminskyqPLMScript} creates an array-format qPLMarr object using
#' the standard output images from
#' \href{http://cad4.cpac.washington.edu/ROTOPOLhome/ROTOPOL.htm}{Rotopol}(tm)
#' software. Images are selected from the working
#' directory using \code{Sys.glob}, allowing this function to build qPLMarr
#' objects without user interaction (e.g, in an unattended "overnight" script).
#'
#' @details
#' More information on the Rotopol(tm) system is available at
#' \url{http://cad4.cpac.washington.edu/ROTOPOLhome/ROTOPOL.htm}.
#'
#' @param sample.name A name that will be used as a stem to label output
#' files from this function (images, etc.).
#'
#' @param bitmap.glob string to direct Sys.glob() to I, |sind|, and phi bitmaps.
#'
#' @param mask.glob string to direct Sys.glob() to mask bitmap.
#'
#' @param north.thickness Numerical; the thickness of the (masked) specimen at
#' its "North" (top center) edge. If no mask is specified, this thickness
#' should be measured at the top center edge of the image.
#'
#' N, S, W, & E thicknesses can be specified separately to model 'wedged'
#' specimens. If the specimen can be assumed to have a constant thickness
#' across the field of view, the N value can be measured anywhere on the
#' image, and will be used to fill in S, W, & E by default.
#'
#' @param south.thickness,west.thickness,east.thickness Numerical; optional
#' thickness measurements from each corner of the
#'
#' @param wavelength Numerical; the center of the illumination wavelength.
#' Default is a 532nm green filter.
#'
#' @param birefringence Numerical; the expected birefringence of the specimen.
#' Default value is an empirical lab standard of 0.0005 for typical bone.
#' Images with several materials of varying birefringence currently
#' require multiple input runs with separate masks for each material.
#'
#' @param pixel Numerical; the absolute size of a pixel in the specimen plane.
#' Intended to be in microns, but the value is arbitrary and not
#' currently linked to specific dimensions in other functions. Default
#' value is lab standard for a Leica A6 macroscope at 1.25x.
#'
#' @param up Numerical; a 'map-view' compass orientation to track a specific axis.
#' Can be used, for example, to specify the dorsal side of a long bone
#' cross-section. This value is passed to the qPLM object's attributes
#' and does not influence the values in the result array.
#'
#' @param mask Boolean; should a mask image be used to define a region of
#' interest? Mask images must be binary (black and white), and have the
#' same pixel dimensions as the other input images. Black pixels (value of
#' 0) will be dropped from the output object, and the output array will
#' be cropped down to a rectangular 10-pixel border around the ROI.
#'
#' @param data.type Character string; description of section type in broad
#' terms, e.g., <<"diaphyseal cross section">>. Keeps specimen data with
#' qPLMarr object as an attribute, and can also be used for watchdogs in
#' analysis functions. For example, applying \code{centroidCorr} to an object
#' will add <<"centroid-corrected">> to the data.type attribute. In future
#' versions, \code{centroidCorr} may prompt the user for confirmation if
#' it is applied to an object that is already labeled as centroid-corrected.
#' Default value is "generic."
#'
#' @return A qPLMarr object that records pixel-by-pixel estimates of I
#' (transmissivity), phi (angular orientation in plane), and theta
#' (angular orientation out of plane). qPLMarr objects also keep
#' relevant data regarding the specimen: thickness, illumination
#' wavelength, the birefringence parameter used to calculate
#' orientations, pixel scale, etc.
#'
#'
# potential tweaks: move "pics.out" routine into a separate function.

KaminskyqPLMScript<-function(sample.name,
                       bitmap.glob,
                       mask.glob=NULL,
                       north.thickness,
                       south.thickness,
                       west.thickness,
                       east.thickness,
                       wavelength=532,
                       birefringence=0.005,
                       pixel=7.5832259,
                       up=0,
                       mask=TRUE,
                       data.type="generic") {
  #require(bmp)
  #requireNamespace(colorspace)
  #require(EBImage)

  Kaminsky.raw<-vector("list", 7)
  # create an object to hold Kaminsky bitmap values
  # and other measured/estimated parameters

  bmpfiles<-Sys.glob(bitmap.glob)
  # Kaminsky bitmap file selection

  if (mask) {
    maskfile<-Sys.glob(mask.glob)
  }
  # mask bitmap file selection

  bmpfiles<-bmpfiles[order(bmpfiles)]
  # slots Kaminsky a0-a2 bitmaps into numerical order

  print("KaminskyqPLM: reading I")
  Kaminsky.raw[[1]]<-bmp::read.bmp(bmpfiles[1])
  # read a0 bitmap (transmittance, I)

  print("KaminskyqPLM: reading |sin d|")
  Kaminsky.raw[[2]]<-bmp::read.bmp(bmpfiles[2])
  # read a1 bitmap (retardance, |sin d|)

  print("KaminskyqPLM: reading phi")
  Kaminsky.raw[[3]]<-bmp::read.bmp(bmpfiles[3])
  # read a2 bitmap (azimuth, phi)

  print("KaminskyqPLM: applying mask")

  ifelse(mask,
         Kaminsky.raw[[7]]<-bmp::read.bmp(maskfile)/255,
         Kaminsky.raw[[7]]<-array(1,dim=dim(Kaminsky.raw[[1]])))
  # read mask or create blank mask object

  dindx<-which(Kaminsky.raw[[7]][,,2]>1/255, arr.ind=TRUE)
  # index of masked pixels

  ifelse(as.logical(min(dindx[,1])<10),
         umin<-1,
         umin<-min(dindx[,1])-10)
  ifelse(as.logical(max(dindx[,1])>dim(Kaminsky.raw[[7]])[1]-10),
         umax<-dim(Kaminsky.raw[[7]])[1],
         umax<-max(dindx[,1])+10)
  ifelse(as.logical(min(dindx[,2])<10),
         vmin<-1,
         vmin<-min(dindx[,2])-10)
  ifelse(as.logical(max(dindx[,2])>dim(Kaminsky.raw[[7]])[2]-10),
         vmax<-dim(Kaminsky.raw[[7]])[2],
         vmax<-max(dindx[,2])+10)

  Kaminsky.raw[[7]]<-as.matrix(Kaminsky.raw[[7]][umin:umax, vmin:vmax,2])
  Kaminsky.raw[[1]]<-as.matrix(Kaminsky.raw[[1]][umin:umax, vmin:vmax,2]
                              *Kaminsky.raw[[7]])
  Kaminsky.raw[[2]]<-as.matrix(Kaminsky.raw[[2]][umin:umax, vmin:vmax,2]
                              *Kaminsky.raw[[7]])
  Kaminsky.raw[[3]]<-as.matrix(Kaminsky.raw[[3]][umin:umax, vmin:vmax,2]
                              *Kaminsky.raw[[7]])
  # trim image arrays to non-zero pixels with a 10-pixel blank border

  Kaminsky.distilled<-array(data=NA,
                           c(dim(Kaminsky.raw[[1]])[2],
                             dim(Kaminsky.raw[[1]])[1], 3))
  # create output object array

  print("KaminskyqPLM: applying parameters")

  thickness<-as.numeric(c(north.thickness,
                          south.thickness,
                          west.thickness,
                          east.thickness))
  Kaminsky.raw[[5]]<-as.numeric(wavelength)
  Kaminsky.raw[[6]]<-as.numeric(birefringence)
  # parameters from arguments

  mid.u.pos<-ceiling((umax-umin+1)/2)
  mid.v.pos<-ceiling((vmax-vmin+1)/2)
  # center N, S, W, E positions for trimmed arrays

  wedge<-matrix(c(1, mid.v.pos,
                  umax-umin+1, mid.v.pos,
                  mid.u.pos, 1,
                  mid.u.pos, vmax-vmin+1),
                byrow=TRUE, nrow=4,
                dimnames=list(c("N","S", "W", "E"), c("u", "v")))
  # N, S, W, E positions in (u, v) coordinates

  wedge<-as.data.frame(cbind(thickness, wedge))
  # data frame for determining wedging

  wedge.f<-lm(thickness ~ u+v, data=wedge)
  # linear fit for wedging

  u.v.pos<-expand.grid(1:(umax-umin+1), 1:(vmax-vmin+1))
  colnames(u.v.pos)<-c("u", "v")
  # pixel positions for by-pixel thickness estimate

  Kaminsky.raw[[4]]<-matrix(predict(wedge.f, u.v.pos), nrow=umax-umin+1)
  # by-pixel thickness estimates

  print("KaminskyqPLM: writing array: I")
  Kaminsky.distilled[,,1]<-t(Kaminsky.raw[[1]]/255)
  # transmittance pixels scaled to 0-1 range

  print("KaminskyqPLM: writing array: Theta")
  Kaminsky.distilled[,,2]<-t(asin
                             (sqrt
                                  (
                                    (asin(Kaminsky.raw[[2]]/255))
                                    *Kaminsky.raw[[5]]
                                    /(2*pi*Kaminsky.raw[[4]]
                                      *1000
                                      *Kaminsky.raw[[6]])
                                    )))
  # retardance pixels transformed to elevation angle (theta)
  # mapped linearly to 0-1 range

  print("KaminskyqPLM: writing array: Phi")
  Kaminsky.distilled[,,3]<-t(Kaminsky.raw[[3]]/255)
  # azimuth pixels (0 degrees to 179 degrees) scaled to 0-1 range


  attributes(Kaminsky.distilled)<-list(dim=dim(Kaminsky.distilled),
                                      dimnames=list(paste("u",seq(length.out=length(Kaminsky.distilled[,1,1])), sep=""),
                                                    paste("v", seq(length.out=length(Kaminsky.distilled[1,,1])), sep=""),
                                                    c("Trans","Theta","Phi")),
                                      thickness_um=Kaminsky.raw[[4]],
                                      wavelength_nm=wavelength,
                                      birefringence=birefringence,
                                      pixel.size_um=pixel,
                                      ccw.skew_deg=up,
                                      dtype=data.type,
                                      class="qPLMarr")
  # sets attributes for qPLMarr object

  rm(Kaminsky.raw)
  gc()

  print("KaminskyqPLM: done!")

  invisible(Kaminsky.distilled)
  return(Kaminsky.distilled)
  # sends calibrated qPLM data to specified object

  gc()
  # may be overkill, but see above
}


