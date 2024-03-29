% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/densBiref.R
\name{densBiref}
\alias{densBiref}
\title{Estimate Birefringence of Bone Tissue From Mineral Density}
\usage{
densBiref(density, dHA = 3.16, birefHA = -0.004, birefC = 0.008933)
}
\arguments{
\item{density}{estimate of mineral density in the sample, in units of
mgHA/ccm.}

\item{dHA}{density of hydroxyapatite in SI units. Default is 3.16.}

\item{birefHA}{birefringence of hydroxyapatite. Default is -4e-3.}

\item{birefC}{birefringence of collagen. Default is 8.933e-3.}
}
\value{
Numeric value giving the estimated birefringence.
}
\description{
\code{densBiref} uses a simple analytical approach to estimating
  birefringence from bone mineral density, as expressed in "milligrams of
  hydroxyapatite per cubic centimeter" (mgHA/ccm).
}
\details{
The actual relationship between bone mineral density and birefringence
  is complicated, because much of the birefringence of collagen is form birefringence
  (Dallemagne and Mélon 1946; Ascenzi and Bonucci 1961), but this approach provides a
  useful heuristic.
}
\references{
Ascenzi, A., & Bonucci, E. (1961). A quantitative investigation of the
  birefringence of the osteon. \emph{Acta Anatomica}, 44(3), 236-262.

  Dallemagne, M. J., & Melon, J. (1946). Nouvelles recherches relatives aux propriétés
  optiques de l'os: La biréfringence de l'os minéralisé; Relations entre les fractions
  organiques et inorganiques de l'os. \emph{Journal of the Washington Academy of Sciences},
  36(6), 181-195.
}
\seealso{
Other qPLM Input Functions: 
\code{\link{buildqPLM}()}
}
\concept{qPLM Input Functions}
