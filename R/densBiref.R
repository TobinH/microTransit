#' @title Estimate Birefringence of Bone Tissue From Mineral Density
#'
#' @description \code{densBiref} uses a simple analytical approach to estimating
#'   birefringence from bone mineral density, as expressed in "milligrams of
#'   hydroxyapatite per cubic centimeter" (mgHA/ccm).
#'
#' @details The actual relationship between bone mineral density and birefringence
#'   is complicated, because much of the birefringence of collagen is form birefringence
#'   (Dallemagne and Mélon 1946; Ascenzi and Bonucci 1961), but this approach provides a
#'   useful heuristic.
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
#' @return Numeric value giving the estimated birefringence.
#'
#' @references Ascenzi, A., & Bonucci, E. (1961). A quantitative investigation of the
#'   birefringence of the osteon. \emph{Acta Anatomica}, 44(3), 236-262.
#'
#'   Dallemagne, M. J., & Melon, J. (1946). Nouvelles recherches relatives aux propriétés
#'   optiques de l'os: La biréfringence de l'os minéralisé; Relations entre les fractions
#'   organiques et inorganiques de l'os. \emph{Journal of the Washington Academy of Sciences},
#'   36(6), 181-195.
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
