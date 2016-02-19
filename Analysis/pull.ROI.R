# subsampling qPLM objects
# arguments:
#  qPLM: a qPLM object
#  layer: value from 1-3, 1 = I, 2 = theta, 3 = phi
pull.ROI<-function(qPLM, layer){
  require(EBImage)
  sel<-"n"
  while (sel != "y") {
    display(qPLM[,,layer], method="raster")
    roi<-locator(2)
    rect(roi$x[1], roi$y[2], roi$x[2], roi$y[1], border="red")
    sub<-qPLM[roi$x[1]:roi$x[2],roi$y[1]:roi$y[2],]
    display(sub[,,layer], method="raster")
    sel<-readline("Keep ROI (y/n)? > ")
  }
  return(sub)
}