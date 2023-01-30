#' @title Multilabel approach for the Tomek Link undersampling algorithm (MLTL)
#'
#' @description This function implements the MLTL algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets, whose aim is to identify
#' tomek links (majoritary instances with a very different neighbor), and remove them. It's like MLeNN, with the number of neighbors being 1.
#'
#' @source Rodolfo M. Pereira, Yandre M.G. Costa, and Carlos N. Silla Jr. MLTL: A multi-label approach for the Tomek Link undersampling algorithm. Neurocomputing, 383:95–105, 2020. ISSN 0925-2312. doi:https://doi.org/10.1016/j.neucom.2019.11.076
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param TH threshold for the Hamming Distance in order to consider an instance different to another one.
#'
#' @return An mldr object containing the preprocessed multilabel dataset
#' @examples
#' \dontrun{
#' library(mldr)
#' MLTL(bibtex, 0.5)
#' }
#' @export
MLTL <- function(D, TH) {

  mldr.resampling::MLeNN(D, TH, 1)

}
