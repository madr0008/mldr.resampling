#' @title Randomly deletes instances with majoritary labels
#'
#' @description This function implements the ML-RUS algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
#' whose aim is to identify instances with majoritary labels, and randomly delete them from the original dataset.
#'
#' @source Charte, F., Rivera, A. J., del Jesus, M. J., & Herrera, F. (2015). Addressing imbalance in multilabel classification: Measures and random resampling algorithms. Neurocomputing, 163, 3-16.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param P Percentage in which the original dataset is increased
#'
#' @return A mld object containing the preprocessed multilabel dataset
#' @examples
#' library(mldr)
#' MLRUS(birds, 25)
#'
#' @export
MLRUS <- function(D, P) {

  mldr::mldr_from_dataframe(D$dataset[-sample(c(1:D$measures$num.instances)[-as.numeric(rownames(D$dataset[rowSums(D$dataset[D$labels[D$labels$IRLbl > D$measures$meanIR,]$index]) > 0,]))],(D$measures$num.instances  / 100) * P,replace=TRUE),], D$labels$index, D$attributes, D$name)

}
