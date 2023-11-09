#' @title Model an Embedded Tissue Block With Registration Marks
#'
#' @description \code{makeBlock} creates a "histoBlock" object to register a
#'   physical cross-section (typically a histology slide) within a MicroCT volume.
#'
#' @details This function works together with \code{registerSlide} to provide
#'   the information necessary to position a section back within the volume of its
#'   parent block, e.g. for importing histological images into a microCT volume.
#'
#' @param B Block width on face plane (match to voxel units, usually millimeters).
#'
#' @param H Block height on face plane (match to voxel units, usually millimeters).
#'
#' @param L1,L2 Upper and lower longitudinal notch offsets. For each, this is
#'   the distance from the nearest lateral face to the far side of the
#'   longitudinal notch, as measured on the front or back faces of the block.
#'
#' @param S1,S2 Upper and lower oblique notch offsets. Measured from the nearest
#'   lateral face to the *near* side of the oblique notch. Lx & Sx should be
#'   measured on the same face (front or back).
#'
#' @param W1,W2 Upper and lower oblique notch widths, measured orthogonal to
#'   oblique notch orientation (essentially the width of the notches as they
#'   were machined).
#'
#' @param Theta1,Theta2 Angle between longitudinal and oblique notches upper
#'   face and lower face. Default is 45 degrees.
#'
#' @param point1,point2 Do the upper/lower notches form an arrowhead pointing
#'   toward the front face of the block?
#'
#' @return A object of type "histoBlock."
#'
#' @family Section Registration Functions
#'
#' @export


makeBlock<-function(B, H, # block width and height
                     L1, S1, W1, Theta1=45, point1=TRUE, # upper notch offsets, thickness, angle, point facing front
                     L2, S2, W2, Theta2=-45, point2=TRUE # lower notch offsets, thickness, angle, point facing front
                      ) {
  block<-vector(mode = "list", length=3)
  names(block)<-c("block dimensions",
                  "top face registration marks",
                  "bottom face registration marks")
  block[[1]]<-as.list(c(B, H))
  block[[2]]<-as.list(c(L1, S1, W1, Theta1, as.logical(point1)))
  block[[3]]<-as.list(c(L2, S2, W2, Theta2, as.logical(point2)))
  names(block[[1]])<-c("block width on section face",
                       "block height on section face")
  names(block[[2]])<-c("top longitudinal notch offset",
                       "top oblique notch offset",
                       "top oblique notch width",
                       "top oblique notch angle", # in degrees; face is zero, longitudinal is 90; top view, counterclockwise
                       "top point towards section face" # TRUE if the oblique and longitudinal notches come to a point on the section face side, FALSE if they come to a point on the back of the block)
                        )
  names(block[[3]])<-c("bottom longitudinal notch offset",
                       "bottom oblique notch offset",
                       "bottom oblique notch width",
                       "bottom oblique notch angle", # in degrees; face is zero, longitudinal is 90; TOP view, counterclockwise
                       "bottom point towards section face" # TRUE if the oblique and longitudinal notches come to a point on the section face side, FALSE if they come to a point on the back of the block)
  )
  attr(block, class)<-"histoBlock"
  return(block)
}
