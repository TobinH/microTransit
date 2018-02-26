#' @title Spatial Clustering of Axial Orientation Data from qPLM
#'
#' @description \code{qPLMClust} produces a hierarchical set of groups that
#'   describe contiguous regions of a specimen with slow axis orientations that
#'   are more similar to each other than they are to neighboring regions.
#'
#' @details \code{qPLMClust} uses iterative agglomerative clustering (Lance &
#'   Williams 1966; 1967) based on Euclidean distances between parameters of the
#'   angular central Gaussian distribution (Tyler 1987), constrained by spatial
#'   neighborhood (Legendre & Legendre 2012:839 - 844).
#'
#'   This function provides an automated means of selecting areas of tissue that
#'   have a distinct "fabric" (e.g., mid-cortical vs. endosteal bone in long
#'   bone cross-section). It was specifically written to avoid the multiple
#'   arbitrary judgement calls that must be made when trying to sub-sample a
#'   bone cross-section using several small ROIs to stand in for a "class" of
#'   tissue.
#'
#'   The comparison is blind to the biological significance of orientation
#'   differences--it is up to the user to select an appropriate hierarchical
#'   level of comparison.
#'
#' @param qPLM A \code{qPLMarr} or \code{qPLMtab} object.
#'
#' @param grainSize Starting subsample dimensions in pixels. \code{qPLM.clust}
#'   will cluster square blocks of \code{grainSize} x \code{grainSize} pixels
#'   based on the angular central Gaussian distribution of their orientations.
#'   Default is a 20 by 20 pixel area (400 pixels), which allows a balance
#'   between discriminating smaller tissue features and providing enough
#'   individual pixels for a tight estimate of distribution.
#'
#' @param cutoff Minimum number of measured pixels present for a block to be
#'   processed. If a \code{grainSize} x \code{grainSize} block has fewer than
#'   \code{cutoff} measured pixels, it will be dropped. Intended to guard
#'   against odd block parameters from edge effects.
#'
#' @param criterion default is \code{"theta,"} which defines clusters based on
#'   the magnitude of the minor axis eigenvalue (low - high values reveal more
#'   or less anisotropy, respectively) and the out-of-plane angles of the major
#'   and second axes of the angular central Gaussian distribution of collagen
#'   fiber orientation per \code{grainSize} x \code{grainSize}) block. This is
#'   roughly equivalent to clustering by brightness in circularly polarized
#'   light microscopy, with some added nuance.
#'
#'   Other available criteria:
#'
#'   \code{"anisotropy"}, which defines clusters based solely on the eigenvalue
#'   magnitude of all three axes, and discards information on the orientation of
#'   the axes (experimental, may be useful for defining kinked/bent tendons or
#'   ligaments vs. surrounding tissue);
#'
#'   \code{"custom"}, which allows users to select from the nine independent
#'   angular central Gaussian descriptors: \describe{ \item{[1:3]}{1st - 3rd
#'   axis eigenvalues;} \item{[4:6]}{1st axis Euclidean coordinates (x,y,z);}
#'   \item{[7:9]}{2nd axis Euclidean coordinates (x,y,z).} } To specify which
#'   descriptors to use, the user must create a list object named "customCrit"
#'   that contains the descriptor numbers as listed above. \code{qPLMClust} will
#'   search the local environment for \code{customCrit}.
#'
#' @param thresHold \code{qPLMClust} will plot all groupings that have a
#'   similarity score of less  than this variable. Useful for getting back to
#'   major tissue differences when there are artifacts in the image that drown
#'   out biologically relevant tissue dissimilarity.
#'
#' @param multi Experimental. Estimating the angular central Gaussian
#'   distribution in parallel using \code{foreach} is much faster,
#'   but sometimes the array binding step using \code{abind} hits a length
#'   mismatch and fails. Until the source of the error can be pinned down,
#'   setting \code{multi} to \code{FALSE} reverts to a nested \code{for}
#'   loop to work around the problem.
#'
#' @return Returns a qPLMclust object, which is a list with four components:
#'   \enumerate{ \item $qPLMtab: the pixel-by-pixel orientatin data in
#'   \code{qPLMtab} format. \item $GaussRaster: a \code{RasterBrick}
#'   representation of the angular central Gaussian distribution parameters per
#'   block. \item $groupMatrix: a matrix with rows that represent processed
#'   pixel blocks, and columns that represent successively more inclusive
#'   agglomerative spatial clusters. The cells in each column give the group
#'   membership for each block at that level, and the contents of the bottom
#'   cell in each column give the similarity value for that agglomerative step.
#'   \item $groupRaster: a \code{RasterBrick} representation of the clusters
#'   that merge at lower similarity values than the value specified by the
#'   \code{thresHold} variable. \code{pullCluster} calls on this component to
#'   let the user page through the high-level clusters and select a level that
#'   reflects their specific biological question (e.g., clusters that separate
#'   an area of tendon insertion relative to surrounding periosteal bone). }
#'
#' @references Lance, G.N., and Williams, W.T., 1966. A generalized sorting
#'   strategy for computer classifications. \emph{Nature 212}:218.
#'
#'   Lance, G.N., and Williams, W.T., 1967. A general theory of classificatory
#'   sorting strategies. I. Hierachical systems. \emph{Computer J 9}:373 - 380.
#'
#'   Legendre, P., and Legendre, L., 2012. \emph{Numerical Ecology}. Elsevier.
#'
#'   Tyler, D.E., 1987. Statistical analysis for the angular central Gaussian
#'   distribution on the sphere. \emph{Biometrika, 74(3)}: 579 - 589.
#'
#' @examples
#' #oldwd<-getwd()
#' #setwd(system.file("extdata", package = "microTransit"))
#' #load("testqPLMarr.R")
#' #load("testqPLMtab.R")
#' #setwd(oldwd)
#'
#' testqPLMclust_arr<-qPLMClust(testqPLMarr)
#' testqPLMclust_tab<-qPLMClust(testqPLMtab)
#'
#'
#' @family qPLM Analysis Functions
#'
#' @export
#'


# Hieronymus Lab 5/8/17


qPLMClust<-function(qPLM,
                    cortical=FALSE,
                    grainSize=20,
                    cutoff=20,
                    criterion="theta",
                    thresHold=0.8,
                    multi = TRUE) {

#   mashButton<-"y"
#   if (attr(qPLM, "class") != "qPLMtab"|"qPLMarr"){
#     mashButton<-readline(prompt = "This doesn't look like a qPLM object. Continue anyway (y/n)? > ")
#     if (mashButton != "y"){
#       break
#     }
#   }

  result<-vector("list", 4)
  names(result)<-c("qPLMtab", "GaussRaster", "groupMatrix", "groupRaster")

  if (attr(qPLM, "class")=="qPLMarr"){
    qPLM<-qPLMTabulate(qPLM) #
  }
  if (cortical){
    xtb<-centroidCorr(qPLM)
  } else {
    xtb<-qPLM
  }
  xtb<-cbind(xtb, xtb[,6]%/%grainSize+1, xtb[,7]%/%grainSize+1) # add grainSize x grainSize addresses
  angU<-max(xtb[,10])
  angV<-max(xtb[,11])

  # abind() is giving me some difficult-to-identify grief in the multi-threaded implementation
  # of ACG estimation per block. Set up flag parameter for single-thread to deal with difficult
  # cases for now.

  debugTime<-Sys.time()

  qPLMGauss<-array(dim=c(angU, angV, 9))

  if (multi){
    wrk<-(parallel::detectCores()*3)%/%4 # count up 3/4 of total CPU cores
    if (is.na(wrk)){ # allows for failure of detectCores()
      wrk<-2
    }

    cl<-parallel::makeCluster(wrk) # set up parallel computation cluster
    doParallel::registerDoParallel(cl)

    qPLMin<-
      foreach::foreach(i=1:angU, .combine=function(...){abind::abind(..., along=3)}) %:%
      foreach::foreach(j=1:angV, .combine=function(...){abind::abind(..., along=2)}, .packages = "microTransit") %dopar% {
        tempGauss<-0
        vec<-vector(mode="numeric", length=9)
        vec[]<-NA
        if (length(which(xtb[,10]==i&xtb[,11]==j))>=cutoff){
          tempGauss<-try(angGaussSumm(xtb[which(xtb[,10]==i&xtb[,11]==j), 3:5], tol=1.5e-06))
          vec[1:3]<-tempGauss$lambda.eigval
          vec[4:6]<-tempGauss$lambda.eigvec[,1]
          vec[7:9]<-tempGauss$lambda.eigvec[,2]
        }
        vec
      }

    parallel::stopCluster(cl)

    for(i in 1:angU){
      for(j in 1:angV){
        for(k in 1:9){
          qPLMGauss[i,j,k]<-qPLMin[k,j,i]
        }
      }
    }

  } else {
    # If foreach() leads to headaches on your machine, this is the
    # single-thread nested for() loop that foreach() attempts to
    # execute in parallel:

    for (i in 1:angU){
      for(j in 1:angV){
        tempGauss<-0
        if (length(which(xtb[,8]==i&xtb[,9]==j))>=cutoff){
          tempGauss<-try(ang.Gauss.summ(xtb[which(xtb[,8]==i&xtb[,9]==j), 3:5], tol=1.5e-06))
          qPLMGauss[i,j,1:3]<-tempGauss$lambda.eigval
          qPLMGauss[i,j,4:6]<-tempGauss$lambda.eigvec[,1]
          qPLMGauss[i,j,7:9]<-tempGauss$lambda.eigvec[,2]
          print(paste(i, j, " recorded--", ((i-1)/angU + j/angV/angU)*100, "%"))
        } else {
          print(paste(i, j, " skipped--", ((i-1)/angU + j/angV/angU)*100, "%"))
        }
      }
    }
  }

  debugTime<-Sys.time()-debugTime
  print(debugTime)

  # align all z vectors for first axes as positive
  for (i in 1:angU){
    for (j in 1:angV){
      qPLMGauss[i,j,4:6]<-qPLMGauss[i,j,4:6]*sign(qPLMGauss[i,j,6])
    }
  }

  # align all z vectors for second axes as positive
  for (i in 1:angU){
    for (j in 1:angV){
      qPLMGauss[i,j,7:9]<-qPLMGauss[i,j,7:9]*sign(qPLMGauss[i,j,9])
    }
  }

  gaussBrick<-raster::brick(qPLMGauss)
  names(gaussBrick)<-c("Major_Axis_Eigenvalue",
                       "Second_Axis_Eigenvalue",
                       "Minor_Axis_Eigenvalue",
                       "Major_Axis_x",
                       "Major_Axis_y",
                       "Major_Axis_z",
                       "Second_Axis_x",
                       "Second_Axis_y",
                       "Second_Axis_z")
  raster::plot(gaussBrick) # debug--shows angular central Gaussian distribution parameters per block

  pixCount<-angU*angV # number of blocks


  # spatial similarity matrix
  listW<-spdep::nb2listw(spdep::cell2nb(angV, angU)) # neighbor list
  links.mat.dat<-spdep::listw2mat(listW) # neighbor matrix
  links.mat.dat[which(links.mat.dat>=0.25)]<-1 # make neighbor matrix binary

  # set up pixel dissimilarity matrix
  GaussTab<-matrix(NA, pixCount, 9) # make a container
  rownames(GaussTab)<-c(1:nrow(GaussTab)) # keep track of where the rows belong
  colnames(GaussTab)<-c("Axis1_eval", "Axis2_eval", "Axis3_eval",
                        "Axis1_evec_x", "Axis1_evec_y", "Axis1_evec_z",
                        "Axis2_evec_x", "Axis2_evec_y", "Axis2_evec_z")
  for (i in 1:pixCount){
    GaussTab[i,]<-qPLMGauss[(i-1)%/%angV+1,(i-1)%%angV+1,]
  } # angular central Gaussian distribution per block

  # cut block info down
  if (criterion == "theta"){
    GaussTab<-GaussTab[,c(3,6,9)] # third axis spread, first and second axis theta
  }

  if (criterion == "anisotropy"){
    GaussTab<-GaussTab[,1:3] # first through third axis spread, ignores orientation
  }

  if (criterion == "custom"){
    GaussTab<-GaussTab[,customCrit] # user-specified in local list object "customCrit"
  }

  # dissimilarity of reduced data
  Diff<-dist(GaussTab, method="euclidean")

  # transform to pixel *similarity* matrix
  Diff<-1-Diff
  Diff<-as.matrix(Diff)
  diag(Diff)<-NA

  # drop "NA" cells from pixel matrix (row names hold on to
  #  block positions so the image can be reassembled)
  fullDiff<-which(!is.na(Diff), arr.ind=TRUE)
  Diff<-Diff[unique(fullDiff[,2]), unique(fullDiff[,2])]
  diag(Diff)<-0

  # drop "NA" cells from spatial similarity matrix
  links.mat.dat<-links.mat.dat[unique(fullDiff[,2]), unique(fullDiff[,2])]

  # manage memory a bit...
  rm(fullDiff)
  rm(qPLMin)
  gc()

  # set up objects to hold group membership
  groupList<-rownames(Diff) # groups by rowname
  groupM<-as.matrix(groupList, ncol=1) # matrix for agglomerative steps
  rownames(groupM)<-rownames(Diff)
  groupList<-c(groupList, 1) # add dummy starting similarity value of 1
  groupM<-rbind(groupM, 1) # final row for similarity at each aggl
  rownames(groupM)[nrow(Diff)+1]<-"Similarity"
  groupSize<-vector(mode="numeric", length=nrow(Diff)) # number of "grains" per agglomerated group
  groupSize[]<-1 # initialize each block as an agglomerated group with size 1

  clusterTime<-Sys.time()

  # iterative agglomerative clustering (Lance and Williams? 19xx) constrained
  #  by spatial neighborhood
  while (length(unique(groupM[c(1:nrow(groupM)-1),ncol(groupM)]))>2){

    # constrain pixel similarity by spatial neighbors
    spatSim<-Diff*links.mat.dat
    # ignore lower triangle of the similarity matrix
    spatSim[lower.tri(spatSim)]<-0

    #debug<-raster(spatSim)
    #plot(debug) # an easy way to visualize the similarity matrix

    # pick next group to glom
    test2<-which(spatSim == max(spatSim), arr.ind=TRUE)

    # shuffle group memberships
    groupList[which(groupList == groupList[test2[,2]])]<-groupList[test2[,1]] # update membership
    groupList[length(groupList)]<-max(spatSim) # record similarity value
    groupM<-cbind(groupM, groupList) # add to output matrix

    # update contiguity matrix
    links.mat.dat[test2[,1],]<-links.mat.dat[test2[,1],]+links.mat.dat[test2[,2],]
    links.mat.dat[,test2[,1]]<-links.mat.dat[,test2[,1]]+links.mat.dat[,test2[,2]]
    links.mat.dat[which(links.mat.dat>1, arr.ind = TRUE)]<-1 # coerce to binary
    links.mat.dat[test2[,2],]<-0 # scrub glommed group
    links.mat.dat[,test2[,2]]<-0 # ""

    # update similarity matrix
    alpha<-groupSize[test2[,1]]/(groupSize[test2[,1]]+groupSize[test2[,2]])
    beta<-groupSize[test2[,2]]/(groupSize[test2[,1]]+groupSize[test2[,2]])
    Diff[test2[,1],]<-alpha*Diff[test2[,1],]+beta*Diff[test2[,2],]
    Diff[,test2[,1]]<-alpha*Diff[,test2[,1]]+beta*Diff[,test2[,2]]
    Diff[test2[,2],]<-0 # scrub glommed group
    Diff[,test2[,2]]<-0 # ""
    diag(Diff)<-0 # pull overlapping similarity scores off diagonal

    # update agglomerative group sizes
    groupSize[test2[,1]]<-groupSize[test2[,1]]+groupSize[test2[,2]]
    groupSize[test2[,2]]<-0 # scrub glommed group

    # scrolling numbers for the impatient. Look busy, your PI is coming.
    print(c(ncol(groupM), groupM[nrow(groupM), ncol(groupM)]))

  }

  clusterTime<-Sys.time()-clusterTime
  print(clusterTime)

  # save image representations of highest hierarchical levels
  hierGroup<-raster::raster(nrows=angU, ncols=angV, xmn=0, xmx=1, ymn=0, ymx=1) # make a raster
  spatGroups<-raster::stack(hierGroup)
  for (i in ncol(groupM):1){
    if (groupM[nrow(groupM), i]>thresHold){
      break
    }
    for (j in 1:nrow(groupM)-1){
      hierGroup[as.numeric(groupM[j,1])]<-as.numeric(groupM[j,i])
    }
    spatGroups<-raster::addLayer(spatGroups, hierGroup)
    print(groupM[nrow(groupM),i])
    raster::values(hierGroup)<-NA
  }

  names(spatGroups)<-c(1:raster::nlayers(spatGroups))

  plot(spatGroups)

  result$qPLMtab<-xtb # place pixel-scale data with group addresses in results for later analysis
  result$GaussRaster<-gaussBrick #
  result$groupMatrix<-groupM
  result$groupRaster<-spatGroups

  attr(result, "thickness_um")<-attr(xtb, "thickness_um")
  attr(result, "wavelength_nm")<-attr(xtb, "wavelength_nm")
  attr(result, "birefringence")<-attr(xtb, "birefringence")
  attr(result, "pixel.size_um")<-attr(xtb, "pixel.size_um")
  attr(result, "ccw.skew_deg")<-attr(xtb, "ccw.skew_deg")
  attr(result, "dtype")<-attr(xtb, "dtype")
  attr(result, "class")<-"qPLMclust"

  return(result)

}




