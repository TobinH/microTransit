# stress eigenfunction worksheet

pos<-tadorna.data$pixels[,6:7]
pos<-scale(pos,center=TRUE,scale=FALSE)
mF<-t(pos)%*%pos
stress.models<-array(data=0, dim=c(dim(pos)[1],3,12,5))

s.xyz<-function(x){
  pos<-scale(x,center=TRUE,scale=FALSE)
  mF<-t(pos)%*%pos
  Ixx<-mF[1,1]
  Iyy<-mF[2,2]
  J<-sum(diag(mF))
  stresses<-array(data=0,dim=c(nrow(pos),3,4))
  stresses[,1,]<-pos[,2]/J
  stresses[,2,]<-pos[,1]/J
  stresses[,3,1]<-pos[,1]/Ixx
  stresses[,3,3]<--pos[,1]/Ixx
  stresses[,3,2]<-pos[,2]/Iyy
  stresses[,3,4]<--pos[,2]/Iyy
  return(stresses)
}
