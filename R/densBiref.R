#' @title Estimate Birefringence of Bone Tissue From Mineral Density
#'
#' @description \code{densBiref} uses a naive analytical approach to estimating
#' birefringence from bone mineral density, as expressed in "milligrams of
#' hydroxyapatite per cubic centimeter" (mgHA/ccm).
#'
#' @param density estimate of mineral density in the sample, in units of
#'   mgHA/ccm.
#'
#' @param dHA density of hydroxyapatite in SI units. Default is 3.16.
#'
#' @param birefHA birefringence of hydroxyapatite. Default is -4e-3.
#'
#' @param birefC birefringence of collagen. Default is 8.933e-3.
#'
#' @return numeric value giving the estimated birefringence.
#'
#' @family qPLM Input Functions
#'
#' @export

densBiref<-function(density,
                     dHA = 3.16,
                     birefHA = -4e-3,
                     birefC = 8.933e-3) {
 birefr<-((density/(dHA*1000)) * birefHA) + ((1 - (density/(dHA*1000))) * birefC)
 return(birefr)
}
