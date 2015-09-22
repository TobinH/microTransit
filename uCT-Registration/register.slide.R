# slide positioning from registration marks

register.slide<-function(B, H, # block width and height
                         L1, S1, W1, Theta1=45, # upper notch offsets, thickness, angle
                         L2, S2, W2, Theta2=-45, # lower notch offsets, thickness, angle
                         D1, F1, # upper notch spread and upper angled notch width from slide 
                         D2, F2  # lower notch spread and lower angled notch width from slide
                         ) {
  dir1<-sign(Theta1)
  dir2<-sign(Theta2)
  ifelse(F1>F2, Slide.PsiY<-(Theta1/180*pi)-(acos(F1/W1)*dir1), Slide.PsiY<-(Theta2/180*pi)-(acos(F2/W2)*dir2))
  Slide.z1<-(sin(Slide.PsiY)*D1*dir1)+(cos(Slide.Psi)*D1/tan(Theta1/180*pi)*dir1)-((S1-L1)/tan(Theta1/180*pi)*dir1)
  Slide.z2<-(sin(Slide.PsiY)*D2*dir2)+(cos(Slide.Psi)*D2/tan(Theta2/180*pi)*dir2)-((S2-L2)/tan(Theta2/180*pi)*dir2)
  Slide.cz1<-Slide.z1-(tan(Slide.PsiY)*(B/2-L1)*dir1)
  Slide.cz2<-Slide.z2-(tan(Slide.PsiY)*(B/2-L2)*dir2)
  Slide.PsiX<-asin(Slide.cz2-Slide.cz1)/H
  Slide.centroid.z<-(Slide.cz1+Slide.cz2)/2
  Slide.Euler.X<-cos(atan(tan(Slide.PsiX)/tan(Slide.PsiY)))
  Slide.Euler.Y<-sin(atan(tan(Slide.PsiX)/tan(Slide.PsiY)))
  Slide.Euler.theta<-atan(tan(Slide.PsiY)/Slide.Euler.X)/pi*180
  results<-vector("list", length=4)
  names(results)<-c("Z", "Euler.X", "Euler.Y", "Euler.Angle")
  results<-c(Slide.centroid.z, Slide.Euler.X, Slide.Euler.Y, Slide.Euler.theta)
  return(results)
}