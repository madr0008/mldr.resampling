#' @title Multi-label oversampling based on local label imbalance (MLSOL)
#'
#' @description This function implements the MLSOL algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets, which applies oversampling on difficult regions of the instance space, in order to help classifiers distinguish labels.
#'
#' @source Liu, B., Blekas, K., & Tsoumakas, G. (2022). Multi-label sampling based on local label imbalance. Pattern Recognition, 122, 108294.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param P Percentage in which the original dataset is increased
#' @param k Number of neighbors to be considered when computing the neighbors of an instance
#' @param neighbors Structure with all instances and neighbors in the dataset. If it is empty, it will be calculated by the function
#' @param tableVDM Dataframe object containing previous calculations for faster processing. If it is empty, the algorithm will be slower
#'
#' @return A mld object containing the preprocessed multilabel dataset
#'
#' @export
MLSOL <- function(D, P, k, neighbors=NULL, tableVDM=NULL) {

  minoritary <- unlist(mldr.resampling.env$.mldrApplyFun1(D$labels$freq, function(x) { ifelse(x<0.5,1,0) }, mc.cores=mldr.resampling.env$.numCores))

  d <- c(1:D$measures$num.instances)[D$dataset$.labelcount > 0]

  if (is.null(neighbors)) {
    message("Part 1/3: Calculating neighbors structure")
    neighbors <- mldr.resampling.env$.mldrApplyFun1(getAllNeighbors(D, d, tableVDM), function(x) { x[1:k+1] }, mc.cores=mldr.resampling.env$.numCores)
  } else {
    message("Part 1/3: Neighbors were already calculated. That just saved us a lot of time!")
    neighbors <- mldr.resampling.env$.mldrApplyFun1(neighbors, function(x) { x[0:k] }, mc.cores=mldr.resampling.env$.numCores)
  }

  message("Part 2/3: Calculating auxiliary structures")

  C <- getC(D, d, neighbors, k)

  S <- getS(D, d, C, minoritary)

  w <- getW(D, S)

  t <- initTypes(C, neighbors, k, minoritary, D, d)

  genNum <- (D$measures$num.instances/100) * P

  seedInstances <- sample(d, size=genNum, replace=TRUE, prob=w)

  message("Part 3/3: Generating new instances")

  newSamples <- mldr.resampling.env$.mldrApplyFun2(seedInstances, function(i) {
    generateInstanceMLSOL(i, sample(neighbors[[i]], size=1), t, D)
  }, mc.cores=mldr.resampling.env$.numCores)

  mldr::mldr_from_dataframe(rbind(D$dataset[1:D$measures$num.attributes], mldr.resampling.env$.mldrApplyFun1(stats::setNames(as.data.frame(do.call(rbind, newSamples[-1])), names(D$attributes)), unlist, mc.cores=mldr.resampling.env$.numCores)), D$labels$index, D$attributes, D$name)

}
