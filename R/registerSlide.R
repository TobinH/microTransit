#' @title Register the Position of a Slide in a Notched Block
#'
#' @description \code{registerSlide} uses notch measurements on a slide,
#'   together with the notch information from a block object (see
#'   \code{make.block}) to position an arbitrarily cut section within the
#'   reference frame of a machined embedded block. This allows precision
#'   registration of plastic-embedded ground sections to a volume representation
#'   of the block, e.g. a MicroCT data set.
#'
#' @details forthcoming.
#'
#' @param block A histoBlock object created by \code{make.block}.
#'
#' @param D1,D2 Distance on slide from far edge of longitudinal notch to near
#'   edge of oblique notch, for upper/lower face, respectively.
#'
#' @param F1,F2 Distance on slide from near edge to far edge of the oblique
#'   notch.
#'
#' @return A list specifying (1) the z-axis position of the slide's centroid
#'   "Z", and (2-4) an Euler rotation for the slide's orientation relative to
#'   the block's front face, specified as "Euler.X", "Euler.Y", and
#'   "Euler.Angle".
#'
#' @family Section Registration Functions
#'
#' @export


registerSlide<-function(block, # block object from function make.block.R
                         D1, F1, # upper notch spread and upper angled notch width from slide
                         D2, F2  # lower notch spread and lower angled notch width from slide
                         ) {
  B<-as.numeric(block[[1]][1])
  H<-as.numeric(block[[1]][2])
  L1<-as.numeric(block[[2]][1])
  S1<-as.numeric(block[[2]][2])
  W1<-as.numeric(block[[2]][3])
  Theta1<-as.numeric(block[[2]][4])
  point1<-as.logical(block[[2]][5])
  L2<-as.numeric(block[[3]][1])
  S2<-as.numeric(block[[3]][2])
  W2<-as.numeric(block[[3]][3])
  Theta2<-as.numeric(block[[3]][4])
  point2<-as.logical(block[[3]][5])
  dir1<-sign(Theta1)
  dir2<-sign(Theta2)
  ifelse(point1,zdir1<-1, zdir1<--1)
  ifelse(point2,zdir2<-1, zdir2<--1)
  Slide.Psi<-vector("numeric", length=1)
  ifelse(F1>F2, Slide.PsiY<-(Theta1/180*pi)-(acos(W1/F1)*dir1), Slide.PsiY<-(Theta2/180*pi)-(acos(W2/F2)*dir2))
  a1<-sin(Slide.PsiY)*D1*dir1*zdir1
  a2<-sin(Slide.PsiY)*D2*dir2*zdir2
  b1<-cos(Slide.Psi)*D1/tan(Theta1/180*pi)*dir1*zdir1
  b2<-cos(Slide.Psi)*D2/tan(Theta2/180*pi)*dir2*zdir2
  c1<-(S1-L1)/tan(Theta1/180*pi)*dir1*zdir1
  c2<-(S2-L2)/tan(Theta2/180*pi)*dir1*zdir2
  Slide.z1<-a1+b1-c1
  Slide.z2<-a2+b2-c2
  Slide.cz1<-Slide.z1-(tan(Slide.PsiY)*(B/2-L1)*dir1)
  Slide.cz2<-Slide.z2-(tan(Slide.PsiY)*(B/2-L2)*dir2)
  Slide.PsiX<-asin((Slide.cz2-Slide.cz1)/H)
  Slide.centroid.z<-(Slide.cz1+Slide.cz2)/2
  Slide.Euler.X<-sin(atan(tan(Slide.PsiX)/tan(Slide.PsiY)))
  Slide.Euler.Y<-cos(atan(tan(Slide.PsiX)/tan(Slide.PsiY)))
  Slide.Euler.theta<-atan(tan(Slide.PsiY)/Slide.Euler.Y)/pi*180
  results<-vector("list", length=4)
  names(results)<-c("Z", "Euler.X", "Euler.Y", "Euler.Angle")
  results<-c(Slide.centroid.z, Slide.Euler.X, Slide.Euler.Y, Slide.Euler.theta)
  return(results)
}
