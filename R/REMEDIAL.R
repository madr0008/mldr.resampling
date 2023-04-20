#' @title Decouples highly imbalanced labels
#' @description This function implements the REMEDIAL algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
#' whose aim is to decouple frequent and rare classes appearing in the same instance. For doing so, it aggregates new instances to the dataset
#' and edit the labels present in them.
#' @source F. Charte, A. J. Rivera, M. J. del Jesus, F. Herrera. "Resampling Multilabel Datasets by Decoupling Highly Imbalanced Labels". Proc. 2015 International Conference on Hybrid Artificial Intelligent Systems (HAIS 2015), pp. 489-501, Bilbao, Spain, 2015. Implementation from the original \code{mldr} package
#' @param mld \code{mldr} object with the multilabel dataset to preprocess
#' @return An mldr object containing the preprocessed multilabel dataset
#' @examples
#' library(mldr)
#' REMEDIAL(birds)
#' @export
REMEDIAL <- function(mld) decoupleImbalancedLabels(mld, mld$measures$scumble)

decoupleImbalancedLabels <- function(mld, atkLevel) {
  mldbase <- mld[mld$dataset$.SCUMBLE <= atkLevel]
  mldhigh <- mld[mld$dataset$.SCUMBLE > atkLevel]  # Samples with coocurrence of highly imbalanced labels

  # Indexes of minority and majority labels
  minIndexes <- mld$labels[mld$labels$IRLbl > mld$measures$meanIR, "index"]
  majIndexes <- mld$labels[mld$labels$IRLbl <= mld$measures$meanIR, "index"]

  # Duplicate rows affected by coocurrence of highly imbalanced labels
  ninstances <- mldhigh$measures$num.instances
  mldhigh$dataset[(ninstances+1):(ninstances*2), ] <- mldhigh$dataset
  row.names(mldhigh$dataset) <- 1:(ninstances*2)

  # Decouple majority and minority labels
  mldhigh$dataset[1:ninstances, minIndexes] <- 0
  mldhigh$dataset[(ninstances+1):(ninstances*2), majIndexes] <- 0

  mldbase + mldhigh # Join the instances without changes with the filtered ones
}
