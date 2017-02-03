# user-specified area subsampling of qPLM objects around a central point
# arguments:
#  qPLM: a qPLM object
#  layer: value from 1-3; this is the layer to display for selecting ROI (1 = I, 2 = theta, 3 = phi).
#         *subsample from ROI will include all of the 'stacked' layers from the qPLM object.
#  pixels: approximate number of pixels required for sample (rounding for a square region may pull more than specified)
pull.qPLM.sample<-function(qPLM, layer, pixels){
  require(EBImage)
  frameSize<-ceiling(sqrt(pixels)/2)
  sel<-"n"
  while (sel != "y") {
    display(qPLM[,,layer], method="raster")
    roi<-locator(1)
    rect((roi$x)-frameSize, (roi$y)-frameSize, (roi$x)+frameSize, (roi$y)+frameSize, border="red")
    sub<-as.array(qPLM[((roi$x)-frameSize):((roi$x)+frameSize), ((roi$y)-frameSize):((roi$y)+frameSize),])
    attributes(sub)<-c(attributes(sub), attributes(qPLM)[3:7])
    display(sub[,,layer], method="raster")
    sel<-readline("Keep subsample (y/n)? > ")
  }
  return(sub)
}