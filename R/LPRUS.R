#' @title Randomly deletes instances with majoritary labelsets
#'
#' @description This function implements the LP-RUS algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
#' whose aim is to identify instances with majoritary labelsets, and randomly delete them from the original dataset.
#'
#' @source Charte, F., Rivera, A. J., del Jesus, M. J., & Herrera, F. (2015). Addressing imbalance in multilabel classification: Measures and random resampling algorithms. Neurocomputing, 163, 3-16.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param P Percentage in which the original dataset is increased
#'
#' @return A mld object containing the preprocessed multilabel dataset
#' @examples
#' library(mldr)
#' LPRUS(birds, 25)
#'
#' @export
LPRUS <- function(D, P) {

  #Calculate the number of samples to deleted in order to decrease in percentage P
  samplesToDelete <- (D$measures$num.instances  / 100) * P

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

  #Obtain instances with majoritary labelsets
  majBag <- c()
  i <- length(labelSetBag)
  while (length(labelSetBag[[i]]) > meanSize) {
    majBag[length(majBag) + 1] <- i
    i <- i - 1
  }

  if (is.null(majBag)) {

    message("There is no imbalance in terms of labelsets. The output dataset is not changed")
    D

  } else {

    #Calculate number of instances to delete, and remove them (process majBag from smallest to largest)
    meanRed <- round(samplesToDelete/length(majBag)) #Round o ceiling?
    meanSize <- floor(meanSize)
    remainder <- rep(0, length(majBag))
    aux <- c()
    for (i in length(majBag):1) {
      rBag <- min(length(labelSetBag[[majBag[i]]]) - meanSize, meanRed)
      remainder[i] <- meanRed - rBag
      #Delete instances
      labelSetBag[[majBag[i]]] <- labelSetBag[[majBag[i]]][-sample.int(length(labelSetBag[[majBag[i]]]), size=rBag, replace=FALSE)]
      aux <- c(aux, rep(i, length(labelSetBag[[majBag[i]]])))
    }

    #Distribute among bags
    i <- length(remainder)
    while ((remainder[i] != 0) && (i >= 1)) {
      x <- sample(aux, size = round(remainder[i]), replace=FALSE)
      bags <- table(x)
      for (j in 1:length(bags)) {
        labelSetBag[[majBag[as.numeric(names(bags)[j])]]] <- labelSetBag[[majBag[as.numeric(names(bags)[j])]]][-sample.int(length(labelSetBag[[majBag[as.numeric(names(bags)[j])]]]), size=bags[[j]], replace=FALSE)]
      }
      aux <- vecsets::vsetdiff(aux, x)
      aux <- aux[(length(labelSetBag[[majBag[i]]]) + 1):length(aux)]
      i <- i - 1
    }

    mldr::mldr_from_dataframe(D$dataset[unlist(labelSetBag),], D$labels$index, D$attributes, D$name)

  }

}
