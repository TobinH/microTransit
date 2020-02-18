#' @title Graphical Representations of Angular Central Gaussian Distributions from qPLM
#'
#' @description \code{acgEllipse} creates an interactive three-dimensional figure
#'   portraying the distribution of axial orientation data.
#'
#' @details Uses the \code{ellipsoid3d()} function from the \code{shapes3d} demo
#'   in package \code{rgl} to plot the eigenvalues and eigenvectors of the
#'   lambda matrix, describing the major axes of an Angular Central Gaussian
#'   distribution.
#'
#' @param angGaussObj Output from the \code{\link{angGaussSumm}} function.
#'
#' @return Silent. Called for the side-effect of producing an \code{rgl} object.
#'
#' @references Tyler, D.E., 1987. Statistical analysis for the angular central
#'   Gaussian distribution on the sphere. \emph{Biometrika, 74(3)}: 579 - 589.
#'
#' @family qPLM Illustration Functions
#'
#' @export

#
acgEllipse<-function(angGaussObj){
  # uses "ellipsoid3d()" code from "shapes3d" demo in rgl
  ellipsoid3d <- function(rx,
                          ry,
                          rz,
                          n=60,
                          ctr=c(0,0,0),
                          ...) {
    degvec <- seq(0,pi,length=n)
    ecoord2 <- function(p) {
      c(rx*cos(p[1])*sin(p[2]),ry*sin(p[1])*sin(p[2]),rz*cos(p[2])) }
    v <- apply(expand.grid(2*degvec,degvec),1,ecoord2)
    v <- rbind(v,rep(1,ncol(v)))
    e <- expand.grid(1:(n-1),1:n)
    i1 <- apply(e,1,function(z)z[1]+n*(z[2]-1))
    i2 <- i1+1
    i3 <- (i1+n-1) %% n^2 + 1
    i4 <- (i2+n-1) %% n^2 + 1
    i <- rbind(i1,i2,i4,i3)
    return(rgl::qmesh3d(v,i,material=list(...)))
  }
  ellip<-ellipsoid3d(rx = angGaussObj$lambda.eigval[1],
                   ry = angGaussObj$lambda.eigval[2],
                   rz = angGaussObj$lambda.eigval[3],
                   qmesh = TRUE,
                   n = 60,
                   color = "orange")
  posEigenvectors<-angGaussObj$lambda.eigvec * sign(angGaussObj$lambda.eigvec[3,1]) # positive z axis 1
  rgl::open3d()
  rgl::bg3d("black")
  rgl::axes3d(color = "white")
  rgl::shade3d(rgl::rotate3d(ellip, matrix = posEigenvectors))
  rgl::title3d(main = NULL, xlab = "X", ylab = "Y", zlab = "Z", color = "white")
}
