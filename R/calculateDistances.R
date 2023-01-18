#' Auxiliary function used by to calculate the distances between an instance and the ones with a specific active label
#'
#' @param sample Index of the sample whose distances to other samples we want to know
#' @param rest Indexes of the samples to which we will calculate the distance
#' @param label Label that must be active
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#'
#' @return A list with the distance to the rest of samples
#' @export
#'
#' @examples
#' \dontrun{
#' library(mldr)
#' calculateDistances(25,c(75,65,89,23),300,bibtex)
#' }
calculateDistances <- function(sample, rest, label, D) {

  sapply(rest, function(y) {
    ifelse(y == sample,
           Inf, #In order not to choose its own
           sqrt( #Square root because euclidean distance
             sum( #Summing distances between numeric and non numeric attributes
               sum( #For numeric attributes: square of the difference
                 (D$dataset[sample,D$attributesIndexes[D$attributes[1:D$measures$num.inputs]=="numeric"]] - D$dataset[y,D$attributesIndexes[D$attributes[1:D$measures$num.inputs]=="numeric"]])^2
               ),
               sum( #For non numeric attributes: Value Difference Measure (VDM)
                 sapply(D$attributesIndexes[D$attributes[1:D$measures$num.inputs]!="numeric"], function(x) {
                   table1 <- table((D$dataset[D$dataset[x] == D$dataset[sample,x],])[label])/(table(D$dataset[x])[[D$dataset[sample,x]]])
                   table2 <- table((D$dataset[D$dataset[x] == D$dataset[y,x],])[label])/(table(D$dataset[x])[[D$dataset[y,x]]])
                   sum(
                     abs(ifelse(length(table1 == 1), ifelse(names(table1) == "0", stats::setNames(c(table1, 0), c("0","1")), stats::setNames(c(0, table1), c("0","1"))), table1) - ifelse(length(table2 == 1), ifelse(names(table2) == "0", stats::setNames(c(table2, 0), c("0","1")), stats::setNames(c(0, table2), c("0","1"))), table2))
                   )
                 })
               )
             )
           )
    )
  })

}
