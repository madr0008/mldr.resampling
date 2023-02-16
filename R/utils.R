#' Auxiliary function used to calculate the distances between an instance and the ones with a specific active label
#'
#' @param sample Index of the sample whose distances to other samples we want to know
#' @param rest Indexes of the samples to which we will calculate the distance
#' @param label Label that must be active
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#'
#' @return A list with the distance to the rest of samples
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




#' Auxiliary function used to compute the neighbors of an instance
#'
#' @param sample Index of the sample whose neighbors we want to know
#' @param rest Indexes of the samples among which we will search
#' @param label Label that must be active, in order to calculate the distances
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param k Number of neighbors to be returned
#'
#' @return A vector with the indexes inside rest of the neighbors
#'
#' @examples
#' \dontrun{
#' library(mldr)
#' getNN(25,c(75,65,89,23),300,bibtex,3)
#' }
getNN <- function(sample, rest, label, D, k) {

  order(calculateDistances(sample, rest, label, D))[1:k+1]

}



#' Auxiliary function used by MLSMOTE. Creates a synthetic sample based on values of attributes and labels of its neighbors
#'
#' @param seedInstance Sample we are using as "template"
#' @param refNeigh Reference neighbor
#' @param neighbors Neighbors to take into account
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#'
#' @return A synthetic sample derived from the one passed as a parameter and its neighbors
#' @examples
#' \dontrun{
#' library(mldr)
#' newSample(25,28,c(15,28,62),bibtex)
#' }
newSample <- function(seedInstance, refNeigh, neighbors, D) {

  c(
    lapply(D$attributesIndexes, function(i) { #Attributes
      ifelse(D$attributes[[i]] %in% c("numeric", "Date"),
             D$dataset[seedInstance,i] + (D$dataset[refNeigh,i] - D$dataset[seedInstance,i])*stats::runif(1, 0, 1), #Numeric attributes
             utils::tail(names(sort(table(D$dataset[neighbors, i]))), 1)) #Non numeric attributes
    }),
    sapply(lapply(D$dataset[c(seedInstance, neighbors),D$labels$index], sum), function(x) { #Labels
      ifelse(x > ((length(neighbors)+1)/2), 1, 0)
    }),
    rep(NA, length(D$dataset) - D$measures$num.attributes) #Other measures like labelcount, SCUMBLE
  )

}



