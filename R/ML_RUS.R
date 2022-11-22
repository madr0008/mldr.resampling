#' Title
#'
#' @param D Original multilabel dataset
#' @param P Percentage in which the original dataset is reduced
#'
#' @return Processed multilabel dataset, more balanced
#' @export
#'
#' @examples
#' ML_RUS(bibtex, 25)
ML_RUS <- function(D, P) {

  mldr_from_dataframe(D$dataset[-sample(c(1:D$measures$num.instances)[-as.numeric(rownames(D$dataset[rowSums(D$dataset[D$labels[D$labels$IRLbl > D$measures$meanIR,]$index]) > 0,]))],(D$measures$num.instances  / 100) * P,replace=FALSE),], D$labels$index, D$attributes, D$name)

}
