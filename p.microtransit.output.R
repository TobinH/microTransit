# parallel execution of microtransit output function

p.microtransit.output<-function(x,sample.name) {
  require(parallel)
  require(foreach)
  require(doParallel)
  require(EBImage)
  require(colorspace)
  st.t<-Sys.time()
  trans<-Image(x[,,1],colormode="Grayscale")
  theta<-Image(x[,,2],colormode="Grayscale")
  phi<-Image(x[,,3],colormode="Grayscale")
  writeImage(trans, file=paste(sample.name,"_trans.tif",sep=""), type="tiff")
  writeImage(theta, file=paste(sample.name,"_theta.tif",sep=""), type="tiff")
  writeImage(phi, file=paste(sample.name,"_phi.tif",sep=""), type="tiff")
  percep.coords<-array(data=NA, dim(x))
  percep.coords[,,1]<-imageData(theta)*73.2
  percep.coords[,,2]<-imageData(theta)*57.65
  percep.coords[,,3]<-imageData(phi)*360
  converted.array<-array(data=NA, dim(x))
  inkwell<-polarLUV(0,0,0)
  conv.col<-matrix(data=NA,nrow=dim(x)[1],ncol=3)
  cl<-makeCluster(detectCores(),outfile="")
  #registerDoSEQ()
  registerDoParallel(cl)
  converted.list<-foreach(k=1:dim(x)[2], .inorder=TRUE, .packages="colorspace") %dopar% {
    for (m in 1:dim(x)[1]) {
      inkwell@coords<-as.matrix(percep.coords[m,k,])
      conv.col[m,]<-coords(as(inkwell,"RGB"))
    }
    conv.col
  }
  for (i in  1:dim(x)[2]){
    converted.array[,i,]<-converted.list[[i]]
  }
  converted.image<-Image(converted.array)
  colorMode(converted.image)<-"Color"
  writeImage(converted.image,file=paste(sample.name,"_overview.tif",sep=""), type="tiff")
  stopCluster(cl)
  print(Sys.time()-st.t)
  print("calibrated images written to working directory")
  return()
}
