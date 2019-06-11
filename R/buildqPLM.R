# Hieronymus Lab 6/11/18

#' @title Build qPLMarr Object From Rotating Analyzer Images.
#'
#' @description \code{buildqPLM} creates an array-format qPLMarr object by
#'   combining images of the same static specimen with the rotating analyzer set
#'   at consistent intervals.
#'
#' @details Native R qPLM analysis from multiple images of a static specimen
#'   with a rotating analyzer. This code is an adaptation of the method of
#'   Glazer et al. (1996), with the intent of providing an inexpensive
#'   alternative that can be employed without extensive microscope
#'   customization.
#'
#' @param steps The number of analyzer positions used. Default is six (six
#'   images with 30 degrees rotation of the analyzer between each).
#'
#' @param tiffs.glob A string that identifies any shared part of the file names
#'   in the image series.
#'
#' @param pixelshift \code{buildqPLM} will try to align the input image stack by
#'   translating the images to get the best fit. The algorithm is experimental.
#'   Set this value to \code{0} to skip automated alignment.
#'
#' @param mask.glob A string that identifies the file name of the mask. Required
#'   if \code{mask = TRUE}.
#'
#' @param north.thickness Thickness of the specimen at its 'North' (top) edge.
#'   N, S, W, & E thicknesses can be specified separately to model 'wedged'
#'   specimens that do not have a consistent thickness across the entire
#'   specimen face. If the specimen can be assumed to have a constant thickness
#'   across the field of view, the N value will be used to fill in S, W, & E by
#'   default.
#'
#' @param south.thickness 'South' (bottom) edge thickness.
#'
#' @param west.thickness 'West' (left) edge thickness.
#'
#' @param east.thickness 'East' (right) edge thickness.
#'
#' @param wavelength The center of the illumination wavelength. Default is a
#'   532nm green filter.
#'
#' @param birefringence The expected birefringence of the specimen. Default
#'   value is a semi-empirical lab standard of \code{5e-4} for bone. Images with
#'   several materials of varying birefringence currently require multiple input
#'   runs with separate masks for each material.
#'
#' @param pixel The absolute size of a pixel in the specimen plane. Intended to
#'   be in microns, but the value is arbitrary and not currently linked to
#'   specific dimensions in other functions. Default value is lab standard for
#'   Leica Z6 Macroscope at 1.25x, MA1000 camera with 1x C-mount.
#'
#' @param up 'Map-view' compass orientation to track a specific axis. Can be
#'   used, for example, to specify the dorsal side of a long bone cross-section.
#'   This value is passed to the qPLMarr object's attributes and does not
#'   influence the values in the result array.
#'
#' @param mask If \code{mask = TRUE}, \code{buildqPLM} will use the string in
#'   \code{bitmap.glob} to select a binary mask image. Mask images must have the
#'   same dimensions as the input images. White pixels in the mask correspond to
#'   pixels from the first input image that will be included in the qPLM output;
#'   black pixels are excluded.
#'
#' @param data.type Character string; description of section type in broad
#'   terms, e.g., <<"diaphyseal cross section">>. Keeps specimen data with
#'   qPLMarr object as an attribute, and can also be used for watchdogs in
#'   analysis functions. For example, applying \code{centroidCorr} to an object
#'   will add <<"centroid-corrected">> to the data.type attribute. In future
#'   versions, \code{centroidCorr} may prompt the user for confirmation if it is
#'   applied to an object that is already labeled as centroid-corrected. Default
#'   value is "generic."
#'
#' @return A \code{qPLMarr} object, an three-dimensional array that records
#'   pixel-by-pixel estimates of \strong{I} (transmissivity), \strong{theta}
#'   (angular orientation out of plane), and \strong{phi} (angular orientation
#'   in plane). qPLMarr objects also keep relevant data regarding the specimen:
#'   thickness, illumination wavelength, the birefringence parameter used to
#'   calculate orientations, pixel scale, etc.
#'
#' @references Glazer, A.M., Lewis, J.G., and Kaminsky, W., 1996. An automatic
#'   optical imaging system for birefringent media. \emph{Proc R Soc A
#'   452(1955)}: 2751 - 2765.
#'
#'   Kaminsky, W., Gunn, E., Sours, R., and Kahr, B., 2007. Simultaneous
#'   false-colour imaging of birefringence, extinction and transmittance
#'   at camera speed. \emph{J Microsc 228(2)}: 153 - 164.
#'
#' @examples
#' oldwd<-getwd()
#' setwd(system.file("extdata", package = "microTransit"))
#' testqPLMarr<-buildqPLM(steps = 6, tiffs.glob = "OsteonImg*",
#'                         mask.glob = "OsteonMask*", north.thickness = 40)
#' save(testqPLMarr, "testqPLMarr.R")
#' setwd(oldwd)
#' qPLMTiff("testqPLM", testqPLMarr)
#'
#'
#' @family qPLM Input Functions
#'
#' @export

buildqPLM<-function(steps = 6,
                    tiffs.glob,
                    pixelshift=10,
                    mask.glob=NULL,
                    north.thickness,
                    south.thickness=NULL,
                    west.thickness=NULL,
                    east.thickness=NULL,
                    wavelength=532, # 532nm green filter, lab default
                    birefringence=0.005, # semi-empirical lab standard for bone birefringence
                    pixel=7.5832259, # macroscope default
                    up=0,
                    mask=TRUE,
                    data.type="generic") {

  if(is.null(south.thickness)){
    south.thickness<-north.thickness
    west.thickness<-north.thickness
    east.thickness<-north.thickness
  }

  # debug: timing
  startT<-Sys.time()


  print("Importing files")
  # input file selection by name search
  infiles<-Sys.glob(tiffs.glob)

  # mask file selection by name search
  if (mask) {
    maskfile<-Sys.glob(mask.glob)
  }

  # create an object to hold bitmap values
  # and other measured/estimated parameters
  image.raw<-vector("list", length=7)

  # set tiff input array dimensions
  image.test<-tiff::readTIFF(infiles[1])
  image.raw[[1]]<-array(data=NA,
                        c(dim(image.test)[1],
                          dim(image.test)[2],
                          steps,
                          3))

  # set mask file's index in list once to avoid reshuffling during debugging
  maskIdx<-3

  # set vector of indices for specimen parameters
  parIdx<-c(4:6)

  # set index for array of derived parameters a0 - a2
  derIdx<-7

  # slots input files into numerical order
  infiles<-infiles[order(infiles)]

  # reads input tiffs into list
  for (i in 1:steps){
    image.raw[[1]][,,i,1]<-tiff::readTIFF(infiles[i])[,,1]
  }

  # Experimental: pixel-shift algorithm to deal with refraction error
  #    from rotating polarizers that are not in the specimen plane

  print(Sys.time() - startT)
  placeT<-Sys.time()

  vshiftstart<-pixelshift+1
  vshiftstop<-(dim(image.raw[[1]])[1]-pixelshift)

  ushiftstart<-pixelshift+1
  ushiftstop<-(dim(image.raw[[1]])[2]-pixelshift)


  if (pixelshift != 0){

    print("Aligning Input Images")

    image.raw[[2]]<-array(data=NA,
                          c(dim(image.test)[1]-(2*pixelshift),
                            dim(image.test)[2]-(2*pixelshift),
                            steps,
                            3))

    image.raw[[2]][,,1,]<-image.raw[[1]][vshiftstart:vshiftstop,
                                         ushiftstart:ushiftstop,
                                         1,]

    alignmentArray<-array(data = NA, dim = c(pixelshift*2+1,
                                             pixelshift*2+1,
                                             dim(image.raw[[1]])[3]-1))

    wrk<-(parallel::detectCores()*3)%/%4 # count up 3/4 of total CPU cores
    if (is.na(wrk)){ # allows for failure of detectCores()
      wrk<-2
    }

    cl<-parallel::makeCluster(wrk) # set up parallel computation cluster
    doParallel::registerDoParallel(cl)

    alignmentArray<-
      foreach::foreach(i=1:(dim(image.raw[[1]])[3]-1), .combine=function(...){abind::abind(..., along=3)}) %dopar% {
        alignmentMatrix<-matrix(data=NA, nrow=pixelshift*2+1, ncol=pixelshift*2+1)
        for (j in -pixelshift:pixelshift){
          for (k in -pixelshift:pixelshift){
            alignmentMatrix[(j+pixelshift+1),(k+pixelshift+1)]<-
              sum(((image.raw[[1]][vshiftstart:vshiftstop,
                                   ushiftstart:ushiftstop,i,1])
                   -(image.raw[[1]][(vshiftstart+j):(vshiftstop+j),
                                    (ushiftstart+k):(ushiftstop+k),i+1,1]))^2)
          }
        }
        alignmentMatrix
      }
    # tested single-thread code as backup to foreach() approach:

    # for (i in 1:(dim(image.raw[[1]])[3]-1)){
    #   for (j in -pixelshift:pixelshift){
    #     for (k in -pixelshift:pixelshift){
    #       alignmentArray[(j+pixelshift+1),(k+pixelshift+1),i]<-sum(((image.raw[[1]][vshiftstart:vshiftstop,ushiftstart:ushiftstop,i,1])-(image.raw[[1]][(vshiftstart+j):(vshiftstop+j), (ushiftstart+k):(ushiftstop+k),i+1,1]))^2)
    #     }
    #   }
    # }

    parallel::stopCluster(cl)

    shiftMatrix<-matrix(data=NA, nrow=dim(image.raw[[1]])[3]-1, ncol=2)

    for (m in 1:nrow(shiftMatrix)){
      shiftMatrix[m,]<-which(alignmentArray[,,m]==min(alignmentArray[,,m]), arr.ind=TRUE)-(pixelshift+1)
    }

    for (p in 2:dim(image.raw[[1]])[3]){
      image.raw[[2]][,,p,]<-image.raw[[1]][(vshiftstart+sum(shiftMatrix[1:(p-1),1])):
                                             (vshiftstop+sum(shiftMatrix[1:(p-1),1])),
                                           (ushiftstart+sum(shiftMatrix[1:(p-1),2])):
                                             (ushiftstop+sum(shiftMatrix[1:(p-1),2])),
                                           p,]
    }
  } else {
    image.raw[[2]]<-image.raw[[1]][vshiftstart:vshiftstop,
                                   ushiftstart:ushiftstop,
                                   ,]
  }

  #image.raw[[1]]<-NULL
  #gc()

  # read mask or create blank mask object
  ifelse(mask,
         image.raw[[maskIdx]]<-tiff::readTIFF(maskfile)[vshiftstart:vshiftstop,
                                                        ushiftstart:ushiftstop,
                                                        1], # mask, if present
         image.raw[[maskIdx]]<-array(1,dim=dim(image.raw[[2]][,,1,1]))) # default trimmed frame

  print(Sys.time() - placeT)
  placeT<-Sys.time()
  print("Cropping to mask")

  # index of masked pixels
  dindx<-which(image.raw[[maskIdx]]>1/255, arr.ind=TRUE)

  # set up coordinates for 10-pixel border around masked area
  ifelse(as.logical(min(dindx[,1])<10),
         umin<-1,
         umin<-min(dindx[,1])-10)
  ifelse(as.logical(max(dindx[,1])>dim(image.raw[[maskIdx]])[1]-10),
         umax<-dim(image.raw[[maskIdx]])[1],
         umax<-max(dindx[,1])+10)
  ifelse(as.logical(min(dindx[,2])<10),
         vmin<-1,
         vmin<-min(dindx[,2])-10)
  ifelse(as.logical(max(dindx[,2])>dim(image.raw[[maskIdx]])[2]-10),
         vmax<-dim(image.raw[[maskIdx]])[2],
         vmax<-max(dindx[,2])+10)

  # trim images to 10-pixel border around masked area
  image.raw[[maskIdx]]<-image.raw[[maskIdx]][umin:umax, vmin:vmax]
  image.raw[[2]]<-image.raw[[2]][umin:umax, vmin:vmax,,]

  for (q in 1:steps){
    image.raw[[2]][,,q,1]<-image.raw[[2]][,,q,1]*image.raw[[maskIdx]]
  }

  print(Sys.time() - placeT)
  placeT<-Sys.time()
  print("applying parameters")

  # specimen parameters from arguments
  thickness<-as.numeric(c(north.thickness,
                          south.thickness,
                          west.thickness,
                          east.thickness))
  image.raw[[parIdx[2]]]<-as.numeric(wavelength)
  image.raw[[parIdx[3]]]<-as.numeric(birefringence)

  # center N, S, W, E positions for trimmed arrays
  mid.u.pos<-ceiling((umax-umin+1)/2)
  mid.v.pos<-ceiling((vmax-vmin+1)/2)

  # N, S, W, E positions in (u, v) coordinates
  wedge<-matrix(c(1, mid.v.pos,
                  umax-umin+1, mid.v.pos,
                  mid.u.pos, 1,
                  mid.u.pos, vmax-vmin+1),
                byrow=TRUE, nrow=4,
                dimnames=list(c("N","S", "W", "E"), c("u", "v")))

  # data frame for determining wedging
  wedge<-as.data.frame(cbind(thickness, wedge))

  # linear fit for wedging
  wedge.f<-lm(thickness ~ u+v, data=wedge)

  # pixel positions for by-pixel thickness estimate
  u.v.pos<-expand.grid(1:(umax-umin+1), 1:(vmax-vmin+1))
  colnames(u.v.pos)<-c("u", "v")

  # by-pixel thickness estimates
  image.raw[[parIdx[1]]]<-matrix(predict(wedge.f, u.v.pos), nrow=umax-umin+1)


  print(Sys.time() - placeT)
  placeT<-Sys.time()
  print("a0, a1, a2 setup")

  # stack for each image
  for (i in 1:steps){
    a.phi<-(i-1)*(pi/steps)
    image.raw[[2]][,,i,]<-abind::abind(image.raw[[2]][,,i,1],
                                       image.raw[[2]][,,i,1]*sin(2*a.phi),
                                       image.raw[[2]][,,i,1]*cos(2*a.phi),
                                       along = 3)
  }

  # set up array for parameters a0, a1, a2
  image.raw[[derIdx]]<-array(data=NA,
                             c(dim(image.raw[[2]])[1],
                               dim(image.raw[[2]])[2], 3))

  print(Sys.time() - placeT)
  placeT<-Sys.time()
  print("a0")

  # fill a0: sum of signal over steps
  for (i in 1:length(image.raw[[2]][,1,1,1])){
    for (j in 1:length(image.raw[[2]][1,,1,1])){
      image.raw[[derIdx]][i,j,1]<-sum(image.raw[[2]][i,j,,1])/steps
    }
  }

  print(Sys.time() - placeT)
  placeT<-Sys.time()
  print("a1")

  #fill a1: sum of signal*sin(2*a.phi) over (steps/2)
  for (i in 1:length(image.raw[[2]][,1,1,1])){
    for (j in 1:length(image.raw[[2]][1,,1,1])){
      image.raw[[derIdx]][i,j,2]<-sum(image.raw[[2]][i,j,,2])*(2/steps)
    }
  }

  print(Sys.time() - placeT)
  placeT<-Sys.time()
  print("a2")

  #fill a2: sum of signal*cos(2*a.phi) over (steps/2)
  for (i in 1:length(image.raw[[2]][,1,1,1])){
    for (j in 1:length(image.raw[[2]][1,,1,1])){
      image.raw[[derIdx]][i,j,3]<-sum(image.raw[[2]][i,j,,3])*(2/steps)
    }
  }

  print(Sys.time() - placeT)
  placeT<-Sys.time()
  print("output setup")

  # create output object array
  qPLM.distilled<-array(data=NA,
                        c(dim(image.raw[[2]])[2],
                          dim(image.raw[[2]])[1], 3))

  print(Sys.time() - placeT)
  placeT<-Sys.time()
  print("|sin d|")

  # fill |sin(d)| layer
  for (i in 1:length(qPLM.distilled[1,,1])){
    for (j in 1:length(qPLM.distilled[,1,1])){
      qPLM.distilled[j,i,2]<-(sqrt
                              ((image.raw[[derIdx]][i,j,2]^2)
                                +(image.raw[[derIdx]][i,j,3]^2)))/image.raw[[derIdx]][i,j,1]
    }
  }

  print(Sys.time() - placeT)
  placeT<-Sys.time()
  print("theta")

  # refill as theta layer
  qPLM.distilled[,,2]<-t(asin
                         (sqrt
                           (asin(t(qPLM.distilled[,,2])) # |sin d|
                             *image.raw[[parIdx[2]]] # wavelength
                             /(2*pi*image.raw[[parIdx[1]]] # thickness
                               *1000
                               *image.raw[[parIdx[3]]]) # birefringence
                           )))/(pi/2)

  print(Sys.time() - placeT)
  placeT<-Sys.time()
  print("phi")

  for (i in 1:length(qPLM.distilled[1,,1])){
    for (j in 1:length(qPLM.distilled[,1,1])){
      qPLM.distilled[j,i,3]<-(pi/2) + sign(image.raw[[derIdx]][i,j,3]) * 0.5 * acos(-image.raw[[derIdx]][i,j,2]/sqrt(image.raw[[derIdx]][i,j,2]^2+image.raw[[derIdx]][i,j,3]^2))
      # updated 6/8/18 with criterion from Kaminsky et al 2007 (MilliView description)
    }
  }

  qPLM.distilled[,,3]<-(qPLM.distilled[,,3]/pi)%%1

  print(Sys.time() - placeT)
  placeT<-Sys.time()
  print("I")

  # fill Io layer
  for (i in 1:length(qPLM.distilled[1,,1])){
    for (j in 1:length(qPLM.distilled[,1,1])){
      qPLM.distilled[j,i,1]<-2*image.raw[[derIdx]][i,j,1]
    }
  }

  # sets attributes for qPLMarr object
  attributes(qPLM.distilled)<-list(dim=dim(qPLM.distilled),
                                   dimnames=list(paste("u",seq(length.out=length(qPLM.distilled[,1,1])), sep=""),
                                                 paste("v", seq(length.out=length(qPLM.distilled[1,,1])), sep=""),
                                                 c("Trans","Theta","Phi")),
                                   thickness_um=image.raw[[parIdx[1]]],
                                   wavelength_nm=wavelength,
                                   birefringence=birefringence,
                                   pixel.size_um=pixel,
                                   ccw.skew_deg=up,
                                   dtype=data.type,
                                   class="qPLMarr")
  print(Sys.time() - placeT)
  print(Sys.time() - startT)

  invisible(qPLM.distilled)
  return(qPLM.distilled)
}
