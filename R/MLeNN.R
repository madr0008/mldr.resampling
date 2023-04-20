#' @title Multilabel edited Nearest Neighbor (MLeNN)
#'
#' @description This function implements the MLeNN algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
#' whose aim is to identify instances with majoritary labels, and remove its neihgbors which are too different to them, in terms of active labels.
#'
#' @source Francisco Charte, Antonio J. Rivera, María J. del Jesus, and Francisco Herrera. MLeNN: A First Approach to Heuristic Multilabel Undersampling. Intelligent Data Engineering and Automated Learning -- IDEAL 2014. ISBN 978-3-319-10840-7.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param TH threshold for the Hamming Distance in order to consider an instance different to another one. Defaults to 0.5.
#' @param k number of nearest neighbours to check for each instance. Defaults to 3.
#' @param neighbors Structure with instances and neighbors. If it is empty, it will be calculated by the function
#' @param tableVDM Dataframe object containing previous calculations for faster processing. If it is empty, the algorithm will be slower
#'
#' @return An mldr object containing the preprocessed multilabel dataset
#'
#' @export
MLeNN <- function(D, TH=0.5, k=3, neighbors=NULL, tableVDM=NULL) {

  majBag <- unique(unlist(mldr.resampling.env$.mldrApplyFun1(D$labels[D$labels$IRLbl < D$measures$meanIR,]$index, function(x) { c(1:D$measures$num.instances)[D$dataset[x]==1] }, mc.cores=mldr.resampling.env$mldr.resampling.env$.numCores)))

  toDelete <- ifelse(is.null(neighbors),

                unlist(mldr.resampling.env$.mldrApplyFun2(majBag, function(x) {
                  activeLabels <- D$labels[which(D$dataset[x,D$labels$index] %in% 1),1]
                  neighbors <- getNN(x, c(1:D$measures$num.instances), ifelse(length(activeLabels)==1,activeLabels,sample(activeLabels,1)), D, tableVDM)[1:k+1]
                  numDifferences <- sum(unlist(mldr.resampling.env$.mldrApplyFun1(neighbors, function(y) {
                                                 adjustedHammingDist(x,y,D) > TH
                                               }, mc.cores=mldr.resampling.env$.numCores)))
                  if (numDifferences >= k/2) { x } #Samples to delete
                }, mc.cores=mldr.resampling.env$.numCores))
                ,
                unlist(mldr.resampling.env$.mldrApplyFun2(c(1:length(majBag)), function(x) {
                  n <- neighbors[[x]][1:k]
                  numDifferences <- sum(unlist(mldr.resampling.env$.mldrApplyFun1(n, function(y) {
                    adjustedHammingDist(majBag[[x]],y,D) > TH
                  }, mc.cores=mldr.resampling.env$.numCores)))
                  if (numDifferences >= k/2) { x } #Samples to delete
                }, mc.cores=mldr.resampling.env$.numCores))

              )

  mldr::mldr_from_dataframe(D$dataset[-toDelete[!is.na(toDelete)],], D$labels$index, D$attributes, D$name)

}
