% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qPLMKeymaker.R
\name{qPLMKeymaker}
\alias{qPLMKeymaker}
\title{Create an Orientation Key for Interpreting qPLM Overview Images}
\usage{
qPLMKeymaker(pixel.size)
}
\arguments{
\item{pixel.size}{Dimension in pixels for square output image.}
}
\value{
Silent return in R. Saves an image with the name "x_pixel_key.tif"
  where x is the value of pixel.size, and displays the image to a browser
  window.
}
\description{
\code{qPLMKeymaker} creates a key image for interpreting
combined colatitude (theta) and longitude (phi) "_overview.tif" images.
}
\details{
not yet.
}
\seealso{
Other qPLM Illustration Functions: \code{\link{acgEllipse}},
  \code{\link{qPLMSphereSkinner}}, \code{\link{qPLMTiff}}
}
