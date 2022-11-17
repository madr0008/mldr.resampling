#' Title
#'
#' @param D Original multilabel dataset
#' @param P Percentage in which the original dataset is reduced
#'
#' @return Processed multilabel dataset, more balanced
#' @export
#'
#' @examples
#' LP_RUS(bibtex, 25)
LP_RUS <- function(D, P) {

  #Calculate the number of samples to deleted in order to decrease in percentage P
  samplesToDelete <- (D$measures$num.instances  / 100) * P

  #Create the bag of labelsets, with the samples that have each labelset
  labelSetBag <- rep(list(c()), D$measures$num.labelsets)
  names(labelSetBag) <- names(D$labelset)

  #Obtain the class label of each instance
  classLabels <- mldr_transform(D, type="LP")$classLabel

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

  mldr_from_dataframe(D$dataset[unlist(labelSetBag),], D$labels$index, D$attributes, D$name)

}
