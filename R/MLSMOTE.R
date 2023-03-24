#' @title Synthetic oversampling of multilabel instances (MLSMOTE)
#'
#' @description This function implements the MLSMOTE algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
#' whose aim is to identify instances with minoritary labels, and generate synthetic instances based on their neighbor instances.
#'
#' @source Francisco Charte, Antonio J. Rivera, María J. del Jesus, and Francisco Herrera. MLSMOTE: Approaching imbalanced multilabel learning through synthetic instance generation. Knowledge-Based Systems, 39:385–397, 2015. ISSN 0950-7051. doi:https://doi.org/10.1016/j.knosys.2015.07.019
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param k Number of neighbors to be considered when creating a synthetic instance
#'
#' @return A mld object containing the preprocessed multilabel dataset
#' @examples
#' \dontrun{
#' library(mldr)
#' MLSMOTE(bibtex, 3)
#' }
#' @export
MLSMOTE <- function(D, k) {

  newSamples <- unlist(.mldrApplyFun2(D$labels[D$labels$IRLbl > D$measures$meanIR,]$index, function(x) {
                                  minBag <- c(1:D$measures$num.instances)[D$dataset[x]==1]
                                  .mldrApplyFun1(minBag, function(y) {
                                    neighbors <- minBag[getNN(y, minBag, x, D)[1:k+1]]
                                    refNeigh <- sample(neighbors, size=1)
                                    newSample(y, refNeigh, neighbors, D)
                                  }, mc.cores=.numCores)
                                }, mc.cores=.numCores), recursive = FALSE)

  mldr::mldr_from_dataframe(rbind(D$dataset[1:D$measures$num.attributes], .mldrApplyFun1(stats::setNames(as.data.frame(do.call(rbind, newSamples[-1])), names(D$attributes)), unlist, mc.cores=.numCores)), D$labels$index, D$attributes, D$name)

}
