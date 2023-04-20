#' @title Randomly clones instances with minoritary labelsets
#'
#' @description This function implements the LP-ROS algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
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
#' LPROS(birds, 25)
#'
#' @export
LPROS <- function(D, P) {

  #Calculate the number of samples to generated in order to increase in percentage P
  samplesToIncrease <- (D$measures$num.instances  / 100) * P

  #Create the bag of labelsets, with the samples that have each labelset
  labelSetBag <- rep(list(c()), D$measures$num.labelsets)
  names(labelSetBag) <- names(D$labelset)

  #Obtain the class label of each instance
  classLabels <- mldr::mldr_transform(D, type="LP")$classLabel

  #Include indexes of the instances with a specific labelset in its labelset bag
  for (i in 1:D$measures$num.instances) {
    labelSetBag[[classLabels[i]]][length(labelSetBag[[classLabels[i]]]) + 1] <- i
  }

  #Calculate the mean number of samples per labelset
  meanSize <- D$measures$num.instances/D$measures$num.labelsets

  #Obtain instances with minoritary labelsets
  minBag <- c()
  for (i in 1:length(labelSetBag)) {
    if (length(labelSetBag[[i]]) < meanSize) {
      minBag[length(minBag) + 1] <- i
    } else {
      break
    }
  }

  if (is.null(minBag)) {

    message("There is no imbalance in terms of labelsets. The output dataset is not changed")
    D

  } else {

    #Calculate number of instances to increase, and add them (process minBag from largest to smallest)
    meanInc <- round(samplesToIncrease/length(minBag)) #Round o ceiling?
    remainder <- rep(0, length(minBag))
    aux <- c()
    for (i in length(minBag):1) {
      rBag <- min(meanSize - length(labelSetBag[[minBag[i]]]), meanInc)
      remainder[i] <- meanInc - rBag
      #Clone instances
      labelSetBag[[minBag[i]]] <- c(labelSetBag[[minBag[i]]], labelSetBag[[minBag[i]]][sample.int(length(labelSetBag[[minBag[i]]]), size=rBag, replace=TRUE)])
      aux <- c(aux, rep(i, length(labelSetBag[[minBag[i]]])))
    }

    #Distribute among bags
    remainder[1] <- round(remainder[1] + samplesToIncrease - meanInc*length(minBag))
    i <- 1
    while ((remainder[i] > 0) && (i <= length(remainder))) {
      x <- sample(aux, size = round(remainder[i]), replace=FALSE)
      bags <- table(x)
      for (j in 1:length(bags)) {
        labelSetBag[[minBag[as.numeric(names(bags)[j])]]] <- c(labelSetBag[[minBag[as.numeric(names(bags)[j])]]], labelSetBag[[minBag[as.numeric(names(bags)[j])]]][sample.int(length(labelSetBag[[minBag[as.numeric(names(bags)[j])]]]), size=bags[[j]], replace=FALSE)])
      }
      aux <- vecsets::vsetdiff(aux, x)
      aux <- aux[(length(labelSetBag[[minBag[i]]]) + 1):length(aux)]
      i <- i + 1
    }

    mldr::mldr_from_dataframe(D$dataset[unlist(labelSetBag),], D$labels$index, D$attributes, D$name)

  }

}
