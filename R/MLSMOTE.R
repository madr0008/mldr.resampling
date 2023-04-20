#' @title Synthetic oversampling of multilabel instances (MLSMOTE)
#'
#' @description This function implements the MLSMOTE algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
#' whose aim is to identify instances with minoritary labels, and generate synthetic instances based on their neighbor instances.
#'
#' @source Charte, F., Rivera, A. J., del Jesus, M. J., & Herrera, F. (2015). MLSMOTE: Approaching imbalanced multilabel learning through synthetic instance generation. Knowledge-Based Systems, 89, 385-397.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param k Number of neighbors to be considered when creating a synthetic instance
#' @param tableVDM Dataframe object containing previous calculations for faster processing. If it is empty, the algorithm will be slower
#'
#' @return A mld object containing the preprocessed multilabel dataset
#'
#' @export
MLSMOTE <- function(D, k, tableVDM=NULL) {

  newSamples <- unlist(mldr.resampling.env$.mldrApplyFun2(D$labels[D$labels$IRLbl > D$measures$meanIR,]$index, function(x) {
                                  minBag <- c(1:D$measures$num.instances)[D$dataset[x]==1]
                                  mldr.resampling.env$.mldrApplyFun1(minBag, function(y) {
                                    neighbors <- ifelse(length(minBag) < k, minBag[getNN(y, minBag, x, D, tableVDM)], minBag[getNN(y, minBag, x, D, tableVDM)[1:k+1]])
                                    refNeigh <- sample(neighbors, size=1)
                                    newSample(y, refNeigh, neighbors, D)
                                  }, mc.cores=mldr.resampling.env$.numCores)
                                }, mc.cores=mldr.resampling.env$.numCores), recursive = FALSE)

  mldr::mldr_from_dataframe(rbind(D$dataset[1:D$measures$num.attributes], mldr.resampling.env$.mldrApplyFun1(stats::setNames(as.data.frame(do.call(rbind, newSamples[-1])), names(D$attributes)), unlist, mc.cores=mldr.resampling.env$.numCores)), D$labels$index, D$attributes, D$name)

}
