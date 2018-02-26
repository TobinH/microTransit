#' @title Interactive Selection of Clustered Pixels
#'
#' @description \code{pullCluster} sub-samples a spatially clustered qPLM
#'   dataset (a \code{qPLMclust} object) and returns a \code{qPLMtab} object
#'   with all the pixels in the selected cluster.
#'
#' @details no.
#'
#' @param qPLMclust A \code{qPLMclust} object, returned from \code{qPLMClust}.
#'
#' @examples
#' #oldwd<-getwd()
#' #setwd(system.file("extdata", package = "microTransit"))
#' #load("testqPLMclust.R")
#' #setwd(oldwd)
#'
#' testCluster<-pullCluster(testqPLMclust_arr)
#'
#' @family qPLM Analysis Functions
#'
#' @export


pullCluster<-function(qPLMclust){
  #require(raster)

  # interactive forward/back paging through plot(raster, layer) of
  # result$groupRaster from qPLMClust()

  i<-1
  rl<-" "
  while(rl!="s"){
    raster::plot(qPLMclust$groupRaster, i)
    rl<-readline(prompt = "Select a cluster layer: (f)orward, (b)ack, (s)elect > ")
    if(substr(rl, 1, 1)=="f"){
      i<-i+1
    }
    if(substr(rl, 1, 1)=="b"){
      i<-i-1
    }
    if(i<1){
      i<-1
      print("Can't go back! already at most inclusive clusters.")
    }
  }

  # use raster package's click() to select a cell from plotted values,
  # pass the group value at the appropriate layer

  print("Click cell(s) on the active plot to select cluster(s).")

  print("Selecting the same cluster more than once will *not* result in duplicated pixels, so slips of the mouse can be safely ignored.")

  print("Use the 'finish' item in the graphics window (or the right mouse button context menu) to complete selection.")

  cell<-raster::click(qPLMclust$groupRaster, show=FALSE)[,i]

  cell<-unique(cell) # just in case any clusters were clicked twice by mistake

  # use that value to get a vector of rows via which() on
  # result$groupMatrix at that layer's column

  rowvec<-which(qPLMclust$groupMatrix[,ncol(qPLMclust$groupMatrix)-(i-1)]==cell)

  # use the rows to get a vector of grainSizexgrainSize addresses
  # from result$groupMatrix first column

  addr<-as.numeric(qPLMclust$groupMatrix[rowvec,1])

  # use the addresses to pull rows from result$qPLMtab

  angV<-max(qPLMclust$qPLMtab[,11])

  results<-NULL

  # single-threaded is somewhat slow

  wrk<-(parallel::detectCores()*3)%/%4 # count up 3/4 of total CPU cores
  if (is.na(wrk)){ # allows for failure of detectCores()
    wrk<-2
  }

  cl<-parallel::makeCluster(wrk) # set up parallel computation cluster
  doParallel::registerDoParallel(cl)

  results<-
    foreach::foreach(j=1:length(addr), .combine=rbind) %dopar% {
      U<-(addr[j]-1)%/%angV+1
      V<-(addr[j]-1)%%angV+1
      chunk<-qPLMclust$qPLMtab[which(qPLMclust$qPLMtab[,10] == U & qPLMclust$qPLMtab[,11] == V),]
      chunk
    }

  parallel::stopCluster(cl)

  # single-thread code, just in case:

  # for (i in 1:length(addr)){
  #   U<-(addr[i]-1)%/%angV+1
  #   V<-(addr[i]-1)%%angV+
  #   chunk<-qPLMclust$qPLMtab[which(qPLMclust$qPLMtab[,10] == U & qPLMclust$qPLMtab[,11] == V),]
  #   results<-rbind(results, chunk)
  # }

  # return the resulting matrix
  dimnames(results)<-dimnames(qPLMclust$qPLMtab)
  attr(results, "thickness_um")<-attr(qPLMclust, "thickness_um")
  attr(results, "wavelength_nm")<-attr(qPLMclust, "wavelength_nm")
  attr(results, "birefringence")<-attr(qPLMclust, "birefringence")
  attr(results, "pixel.size_um")<-attr(qPLMclust, "pixel.size_um")
  attr(results, "ccw.skew_deg")<-attr(qPLMclust, "ccw.skew_deg")
  attr(results, "dtype")<-attr(qPLMclust, "dtype")
  attr(results, "class")<-"qPLMtab"

  return(results)

}











