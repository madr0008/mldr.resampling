#' @title Multi-label oversampling based on local label imbalance (MLSOL)
#'
#' @description This function implements the MLSOL algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets, which applies oversampling on difficult regions of the instance space, in order to help classifiers distinguish labels.
#'
#' @source Liu, B., Blekas, K., & Tsoumakas, G. (2022). Multi-label sampling based on local label imbalance. Pattern Recognition, 122, 108294.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param P Percentage in which the original dataset is increased
#' @param k Number of neighbors to be considered when computing the neighbors of an instance
#'
#' @return A mld object containing the preprocessed multilabel dataset
#' @examples
#' \dontrun{
#' library(mldr)
#' MLSOL(bibtex, 3)
#' }
#' @export
MLSOL <- function(D, P, k) {

  minoritary <- unlist(lapply(D$labels$freq, function(x) ifelse(x<0.5,1,0)))

  d <- as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,]))

  neighbors <- getAllNeighbors(D, d, k)

  C <- getC(D, d, neighbors, k)

  S <- getS(D, d, C, minoritary)

  w <- getW(D, d, S)

  t <- initTypes(C, neighbors, k, minoritary)

  genNum <- (D$measures$num.instances/100) * P

  seedInstances <- sample(d, size=genNum, replace=TRUE, prob=w)

  newSamples <- lapply(seedInstances, function(i) {
    generateInstanceMLSOL(i, sample(neighbors[[as.character(i)]], size=1), t, D)
  })

  mldr::mldr_from_dataframe(rbind(D$dataset, lapply(stats::setNames(as.data.frame(do.call(rbind, newSamples[-1])), names(D$dataset)), unlist)), D$labels$index, D$attributes, D$name)

}