# pull analysis data from microtransit object
# note that as of 4 Feb 15 microtransit is not a defined object class
# seems like a natural progression, but one thing at a time

p.qPLMtabulate<-function(x,low.pass=5){
  require(parallel)
  require(foreach)
  require(doParallel)
  st.t<-Sys.time()
  nz.pix<-which(x[,,2]>(low.pass), arr.ind=TRUE)
  # low pass retardance filter to remove "empty" pixels from analysis
  fourx.pos<-matrix(data=NA,nrow=0,ncol=2)
  # empty matrix for 4x4 binned pixel xy positions
  tabulated.data<-matrix(data=NA,nrow=0,ncol=7)
  # empty matrix for [,1:2] angular data from qPLM, [,3:5] xyz transform of same,
  # [,6:7] pixel position in 4x4 bins for later indexing
  # temp.pos<-matrix(nrow=0,ncol=2)
  # temp.ang<-matrix(nrow=0,ncol=7)
  cl<-makeCluster(detectCores(),outfile="")
  # processor count for parallel computing
  registerDoParallel(cl)
  # set up "cluster" of processors
  tab.data<-foreach(i=1:512, combine='rbind') %:%
    foreach(j=1:384, combine='rbind') %dopar% {
      
    }
