#' Auxiliary function used by MLeNN. Computes the Hamming Distance between two instances
#'
#' @param x Index of sample 1
#' @param y Index of sample 2
#' @param D mld \code{mldr} object in which the instances are located
#'
#' @return The Hamming Distance between the instances
#' @examples
#' \dontrun{
#' library(mldr)
#' adjustedHammingDist(1,2,bibtex)
#' }
adjustedHammingDist <- function(x,y,D) {
  length(which(D$dataset[x,D$labels$index] != D$dataset[y,D$labels$index])) / (sum(D$dataset[x,D$labels$index]) + sum(D$dataset[y,D$labels$index]))
}



#' @title Multilabel edited Nearest Neighbor (MLeNN)
#'
#' @description This function implements the MLeNN algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
#' whose aim is to identify instances with minoritary labels, and remove its neihgbors which are too different to them, in terms of active labels.
#'
#' @source Francisco Charte, Antonio J. Rivera, María J. del Jesus, and Francisco Herrera. MLeNN: A First Approach to Heuristic Multilabel Undersampling. Intelligent Data Engineering and Automated Learning -- IDEAL 2014. ISBN 978-3-319-10840-7.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param TH threshold for the Hamming Distance in order to consider an instance different to another one. Defaults to 0.5.
#' @param NN number of nearest neighbours to check for each instance. Defaults to 3.
#'
#' @return An mldr object containing the preprocessed multilabel dataset
#' @examples
#' \dontrun{
#' library(mldr)
#' MLeNN(bibtex, 0.5, 3)
#' }
#' @export
MLeNN <- function(D, TH=0.5, NN=3) {#Obtain indexes of minoritary labels
  minLabels <- D$labels[D$labels$IRLbl > D$measures$meanIR,]$index

  #Obtain indexes of instances with each minority label
  minBag <- unique(unlist(lapply(minLabels, function(x) as.numeric(rownames(D$dataset[D$dataset[x]==1,])))))

  toDelete <- unlist(pbapply::pblapply(minBag, function(x) {
                                                  activeLabels <- D$labels[which(D$dataset[x,D$labels$index] %in% 1),1]
                                                  neighbors <- order(calculateDistances(x, as.numeric(rownames(D$dataset)), ifelse(length(activeLabels)==1,activeLabels,sample(activeLabels,1)), D))[1:NN+1]
                                                  numDifferences <- sum(unlist(lapply(neighbors, function(y) {
                                                    adjustedHammingDist(x,y,D) > TH
                                                  })))
                                                  ifelse(numDifferences >= NN/2,x) #Samples to delete
                                                }))

  mldr::mldr_from_dataframe(D$dataset[-toDelete[!is.na(toDelete)],], D$labels$index, D$attributes, D$name)

}
