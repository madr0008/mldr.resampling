#' @title Randomly clones instances with minoritary labels
#'
#' @description This function implements the ML-ROS algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
#' whose aim is to identify instances with minoritary labels, and randomly clone them.
#'
#' @source Charte, F., Rivera, A. J., del Jesus, M. J., & Herrera, F. (2015). Addressing imbalance in multilabel classification: Measures and random resampling algorithms. Neurocomputing, 163, 3-16.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param P Percentage in which the original dataset is increased
#'
#' @return A mld object containing the preprocessed multilabel dataset
#' @examples
#' library(mldr)
#' library(mldr.resampling)
#' MLROS(birds, 25)
#'
#' @export
MLROS <- function(D, P) {

  #Calculate the number of samples to deleted in order to decrease in percentage P
  samplesToClone <- (D$measures$num.instances  / 100) * P

  #Obtain indexes of minoritary labels
  minLabels <- D$labels[D$labels$IRLbl > D$measures$meanIR,]$index

  #Obtain indexes of instances with each minority label
  minBag <- mldr.resampling.env$.mldrApplyFun1(minLabels, function(x) { as.numeric(rownames(D$dataset[D$dataset[x]==1,])) }, mc.cores=mldr.resampling.env$.numCores)

  #Instance cloning loop
  clonedSamples <- c()
  canClone <- rep(TRUE, length(minBag))
  labelCount <- D$labels$count
  maxCount <- max(D$labels$count)
  while ((samplesToClone > 0) && all(canClone)) {

    samplesToClone <- samplesToClone - sum(canClone)
    clonedSamples <- c(clonedSamples, unlist(mldr.resampling.env$.mldrApplyFun1(c(1:length(minBag)), function(i) {
      if (canClone[i]) {
        return(sample(minBag[[i]],1,replace=FALSE))
        canClone[i] <- !(maxCount/D$labels$count[D$labels$index == minLabels[i]] <= D$measures$meanIR)
      }
    }, mc.cores=mldr.resampling.env$.numCores)))

  }

  mldr::mldr_from_dataframe(D$dataset[c(1:D$measures$num.instances,clonedSamples),], D$labels$index, D$attributes, D$name)

}
