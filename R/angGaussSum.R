#' @title Estimate and Summarize Angular Central Gaussian Distribution for Axial
#'   Data
#'
#' @description \code{angGaussSumm} estimates the angular central Gaussian
#'   distribution of 3D axial data using Tyler's (1987) maximum likelihood
#'   method.
#'
#' @param X An \code{[n,3]} matrix of Euclidean coordinates for axial
#'   orientations, for general use. \code{angGaussSumm} will also parse the
#'   orientations from a \code{qPLMtab} object referenced by this argument,
#'   without requiring specific column addresses.
#'
#' @param alpha Alpha value corresponding to confidence interval \code{1 - CI}.
#'   Default is 95\% CI.
#'
#' @param tol Tolerance for convergence of iterative algorithm. Default value is
#'   the default for the \code{all.equal()} function. Increase this value if the
#'   function returns an identity matrix for the lambda estimate.
#'
#' @return A list with four components: \enumerate{ \item $lambda.matrix: a
#'   matrix of lambda values, \item $lambda.eigval: the eigenvalues of the
#'   lambda matrix, \item $lambda.eigvec: the eigenvectors of the lambda matrix
#'   (major orientations of the ACG distribution), and \item $aCI: angular
#'   confidence intervals for first axis in second axis plane (I-II), first axis
#'   in third axis plane (I-III), and second axis in third axis plane (II-III),
#'   given in both radians and degrees.}
#'
#' @references Tyler, D.E., 1987. Statistical analysis for the angular central
#'   Gaussian distribution on the sphere. \emph{Biometrika, 74(3)}: 579 - 589.
#'
#' @examples
#' #oldwd<-getwd()
#' #setwd(system.file("extdata", package = "microTransit"))
#' #load("testqPLMtab.R")
#' testqPLMGauss_tab<-angGaussSumm(testqPLMtab)
#' save(testqPLMGauss_tab, file = "testqPLMGauss_tab.R")
#' #setwd(oldwd)
#'
#' @family qPLM Analysis Functions
#'
#' @export

angGaussSumm<-function(X, alpha=0.05, tol=sqrt(.Machine$double.eps)){
  if (!is.null(attr(X, "class"))) {
    if (attr(X, "class")=="qPLMtab") {
      X<-X[,3:5]
    }
  }
  lambda<-diag(1, nrow=3, ncol=3) # identity lamba matrix to start
  lambdahat<-matrix(data=0, nrow=3, ncol=3) # lambda estimate for Tyler's eqn (2)
  counter<-0 # number of iterations
  convCrit1<- 3 # container for previous iteration's value of the difference between tr(lambdahat) & q
  convCrit2<- 3 # container for current interations's value of ""
  while (counter < 4 | !isTRUE(all.equal(convCrit1, convCrit2, tolerance = tol))){ # iterative procedure for ML estimate of lambda matrix
    # tolerance can be tweaked for missed convergence in test cases
    # poorly conditioned matrices cause all sorts of problems... check first
    if (rcond(lambda) < tol){
      break
    }
    convCrit1<-convCrit2
    convCrit2<-0
    sumli<-matrix(data=0, nrow=ncol(X), ncol=(ncol(X)))
    sumone<-vector(mode="numeric", length=1)
    sumone<-0
    for (i in 1:nrow(X)){
      sumli<-((X[i,]%*%t(X[i,]))/as.numeric(t(X[i,])%*%solve(lambda)%*%X[i,]))+sumli
      sumone<-(1/as.numeric(t(X[i,])%*%solve(lambda)%*%X[i,]))+sumone
    }
    lambda<-ncol(X)*sumli/sumone
    lambdahat<-ncol(X)/nrow(X)*sumli
    convCrit2<-sum(diag(lambdahat))-ncol(X)
    counter<-counter+1
    cat("Iteration", counter, "*** convergence score =", convCrit2, "*** last convergence score =", convCrit1,
        "*** difference =",  convCrit1-convCrit2, fill=TRUE)

  }
  gam<-eigen(lambda) # distribution axis weights and vectors

  aCI12<-asin(sqrt((nrow(X)*ncol(X))^-1
                         *(ncol(X)+2)^-1
                         *qchisq(1-(alpha/2), df=ncol(X)-1)
                         *gam$values[2]
                         *gam$values[1])
                    /abs(gam$values[2]-gam$values[1])) # confidence interval from Axis I eigenvectors in Axis I-II plane, in radians
  aCI13<-asin(sqrt((nrow(X)*ncol(X))^-1
                         *(ncol(X)+2)^-1
                         *qchisq(1-(alpha/2), df=ncol(X)-1)
                         *gam$values[3]
                         *gam$values[1])
                    /abs(gam$values[3]-gam$values[1])) # confidence interval from Axis I eigenvectors in Axis I-III plane, in radians
  aCI23<-asin(sqrt((nrow(X)*ncol(X))^-1
                         *(ncol(X)+2)^-1
                         *qchisq(1-(alpha/2), df=ncol(X)-1)
                         *gam$values[3]
                         *gam$values[2])
                    /abs(gam$values[3]-gam$values[2])) # confidence interval around Axis II eigenvectors in Axis II-III plane, in radians
  aCI<-matrix(data=0, nrow=2, ncol=3)
  aCI[1,]<-c(aCI12, aCI13, aCI23) # CIs in radians
  aCI[2,]<-aCI[1,]/pi*180 # CIs in degrees
  res<-vector(mode="list", length=4)
  names(res)<-c("lambda.matrix", "lambda.eigval", "lambda.eigvec", "aCI")
  res$lambda.matrix<-lambda
  res$lambda.eigval<-gam$values
  names(res$lambda.eigval)<-c("Axis I", "Axis II", "Axis III")
  res$lambda.eigvec<-gam$vectors
  rownames(res$lambda.eigvec)<-colnames(X)
  colnames(res$lambda.eigvec)<-c("Axis I", "Axis II", "Axis III")
  res$aCI<-aCI
  colnames(res$aCI)<-c("I-II", "I-III", "II-III")
  rownames(res$aCI)<-c("radians", "degrees")
  return(res)
}

