# function to force cortical bone fiber orientations to right-handed (ccw towards viewer)
# or left-handed (cw towards viewer) helices.

# arguments:
#    qPLMtab: a qPLMtab object obtained from qPLMtabulate.R
#    cw: left-handed helix. Default is FALSE.


resolve.theta.ambiguity<-function(qPLMtab,cw=FALSE){
  result<-qPLMtab$pixels
  result[,6:7]<-scale(result[,6:7], scale=FALSE, center=TRUE)
  ifelse (cw, sign<-(-1), sign<-1)
  for (i in 1:nrow(result)){
    if ((result[i,7]*sign)>abs(result[i,6])){
      result[i,3:4]<-result[i,3:4]*-1
    }
    if ((result[i,6]*sign*-1)>abs(result[i,7])){
      if (result[i,1]>(pi/2)){
        result[i,3:4]<-result[i,3:4]*-1
      }
    }
    if ((result[i,6]*sign)>abs(result[i,7])){
      if (result[i,1]<(pi/2)){
        result[i,3:4]<-result[i,3:4]*-1
      }
    }
    if (result[i,4]<0){
      result[i,2]<-result[i,2]+pi
    }
  }
  invisible(result)
  return(result)
}