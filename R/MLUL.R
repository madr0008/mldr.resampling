#' @title Multi-label undersampling based on local label imbalance (MLUL)
#'
#' @description This function implements the MLUL algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets, which applies undersampling, removing difficult instances according to their neighbors.
#'
#' @source Liu, B., Blekas, K., & Tsoumakas, G. (2022). Multi-label sampling based on local label imbalance. Pattern Recognition, 122, 108294.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param P Percentage in which the original dataset is decreased
#' @param k Number of neighbors to be considered when computing the neighbors of an instance
#'
#' @return A mld object containing the preprocessed multilabel dataset
#' @examples
#' \dontrun{
#' library(mldr)
#' MLUL(bibtex, 3)
#' }
#' @export
MLUL <- function(D, P, k) {

  minoritary <- unlist(lapply(D$labels$freq, function(x) ifelse(x<0.5,1,0)))

  d <- as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,]))

  neighbors <- getAllNeighbors(D, d, k)

  rNeighbors <- getAllReverseNeighbors(d, neighbors, k)

  C <- getC(D, d, neighbors, k)

  S <- getS(D, d, C, minoritary)

  w <- getW(D, d, S)

  u <- getU(D, d, rNeighbors, S)

  v <- getV(d, w, u)

  retNum <- (length(d)/100) * (100 - P)

  toKeep <- sample(d, size=ceiling(retNum), replace=FALSE, prob=v)

  mldr::mldr_from_dataframe(D$dataset[toKeep,], D$labels$index, D$attributes, D$name)

}