% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pullCluster.R
\name{pullCluster}
\alias{pullCluster}
\title{Interactive Selection of Clustered Pixels}
\usage{
pullCluster(qPLMclust)
}
\arguments{
\item{qPLMclust}{A \code{qPLMclust} object, returned from \code{qPLMClust}.}
}
\description{
\code{pullCluster} sub-samples a spatially clustered qPLM
  dataset (a \code{qPLMclust} object) and returns a \code{qPLMtab} object
  with all the pixels in the selected cluster node.
}
\details{
This function allows the user to sample an existing spatial clustering
  scheme generated using \code{qPLMClust}. The display starts at the deepest split
  in the cluster analysis, showing the two most dissimilar patches of pixels by
  whichever criteria were chosen for \code{qPLMClust}. Advancing through the nodes by
  entering \code{f} at the prompt shows additional splits in order of increasing
  similarity. Once the desired set of splits is displayed, The user can then enter
  \code{s} at the prompt to select one or more groups of pixels.
}
\examples{
#oldwd<-getwd()
#setwd(system.file("extdata", package = "microTransit"))
#load("testqPLMclust.R")
#setwd(oldwd)

testCluster<-pullCluster(testqPLMclust_arr)

}
\seealso{
Other qPLM Analysis Functions: 
\code{\link{angGaussSumm}()},
\code{\link{centroidCorr}()},
\code{\link{pullROI}()},
\code{\link{pullqPLMSample}()},
\code{\link{qPLMClust}()},
\code{\link{qPLMTabulate}()}
}
\concept{qPLM Analysis Functions}
