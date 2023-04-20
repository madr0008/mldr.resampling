#' @title Multi-label undersampling based on local label imbalance (MLUL)
#'
#' @description This function implements the MLUL algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets, which applies undersampling, removing difficult instances according to their neighbors.
#'
#' @source Liu, B., Blekas, K., & Tsoumakas, G. (2022). Multi-label sampling based on local label imbalance. Pattern Recognition, 122, 108294.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param P Percentage in which the original dataset is decreased
#' @param k Number of neighbors to be considered when computing the neighbors of an instance
#' @param neighbors Structure with all instances and neighbors in the dataset. If it is empty, it will be calculated by the function
#' @param tableVDM Dataframe object containing previous calculations for faster processing. If it is empty, the algorithm will be slower
#'
#' @return A mld object containing the preprocessed multilabel dataset
#'
#' @export
MLUL <- function(D, P, k, neighbors=NULL, tableVDM=NULL) {

  minoritary <- unlist(mldr.resampling.env$.mldrApplyFun1(D$labels$freq, function(x) { ifelse(x<0.5,1,0)} , mc.cores=mldr.resampling.env$.numCores))

  d <- c(1:D$measures$num.instances)[D$dataset$.labelcount > 0]

  if (is.null(neighbors)) {
    message("Part 1/2: Calculating neighbors structure")
    neighbors <- mldr.resampling.env$.mldrApplyFun1(getAllNeighbors(D, d, tableVDM), function(x) { x[1:k+1] }, mc.cores=mldr.resampling.env$.numCores)
  } else {
    message("Part 1/2: Neighbors were already calculated. That just saved us a lot of time!")
    neighbors <- mldr.resampling.env$.mldrApplyFun1(neighbors, function(x) { x[0:k] }, mc.cores=mldr.resampling.env$.numCores)
  }

  message("Part 2/2: Calculating auxiliary structures and removing instances")

  rNeighbors <- getAllReverseNeighbors(d, neighbors, k)

  C <- getC(D, d, neighbors, k)

  S <- getS(D, d, C, minoritary)

  w <- getW(D, S)

  u <- getU(D, d, rNeighbors, S)

  v <- getV(w, u)

  retNum <- (length(d)/100) * (100 - P)

  toKeep <- sample(d, size=ceiling(retNum), replace=FALSE, prob=v)

  mldr::mldr_from_dataframe(D$dataset[toKeep,], D$labels$index, D$attributes, D$name)

}
