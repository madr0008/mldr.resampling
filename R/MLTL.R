#' @title Multilabel approach for the Tomek Link undersampling algorithm (MLTL)
#'
#' @description This function implements the MLTL algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets, whose aim is to identify
#' tomek links (majoritary instances with a very different neighbor), and remove them. It's like MLeNN, with the number of neighbors being 1.
#'
#' @source Pereira, R. M., Costa, Y. M., & Silla Jr, C. N. (2020). MLTL: A multi-label approach for the Tomek Link undersampling algorithm. Neurocomputing, 383, 95-105.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param TH threshold for the Hamming Distance in order to consider an instance different to another one.
#' @param neighbors Structure with instances and neighbors. If it is empty, it will be calculated by the function
#' @param tableVDM Dataframe object containing previous calculations for faster processing. If it is empty, the algorithm will be slower
#'
#' @return An mldr object containing the preprocessed multilabel dataset
#'
#' @export
MLTL <- function(D, TH, neighbors=NULL, tableVDM=NULL) {

  mldr.resampling::MLeNN(D, TH, 1, neighbors, tableVDM)

}
