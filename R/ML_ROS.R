#' Title
#'
#' @param D Original multilabel dataset
#' @param P Percentage in which the original dataset is increased
#'
#' @return Processed multilabel dataset, more balanced
#' @export
#'
#' @examples
#' ML_ROS(bibtex, 25)
ML_ROS <- function(D, P) {

  #Calculate the number of samples to deleted in order to decrease in percentage P
  samplesToClone <- (D$measures$num.instances  / 100) * P

  #Obtain indexes of minoritary labels
  minLabels <- D$labels[D$labels$IRLbl > D$measures$meanIR,]$index

  #Obtain indexes of each minority label
  minBag <- rep(list(c()),length(minLabels))
  #names(minBag) <- minLabels
  for (i in 1:length(minBag)) {
    minBag[[i]] <- as.numeric(rownames(D$dataset[D$dataset[minLabels[i]]==1,]))
  }

  #Instance cloning loop
  clonedSamples <- c()
  canClone <- rep(TRUE, length(minBag))
  labelCount <- D$labels$count
  maxCount <- max(D$labels$count)
  while (samplesToClone > 0) { #&& canClone

    for (i in 1:length(minBag)) {

      if (canClone[i]) {

        clonedSamples[length(clonedSamples) + 1] <- sample(minBag[[i]],1,replace=FALSE)
        samplesToClone <- samplesToClone - 1
        newIRLbl <- maxCount/D$labels$count[D$labels$index == minLabels[i]]
        if (newIRLbl <= D$measures$meanIR) {
          canClone[i] <- FALSE
        }
      }

    }

  }

  mldr_from_dataframe(D$dataset[c(1:D$measures$num.instances,clonedSamples),], D$labels$index, D$attributes, D$name)


}
