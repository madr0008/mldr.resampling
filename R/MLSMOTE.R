#' Auxiliary function used by MLSMOTE. Calculates the distances between an instance and the ones with the same label
#'
#' @param sample Sample whose distances to other samples with the same label we want to know
#' @param minBag Instances with a determined label
#' @param label Minoritary label we want to balance
#'
#' @return A list with the distances to the rest of samples with the same label
#' @export
#'
#' @examples
#' calculateDistances(25,c(75,65,89,23),300)
calculateDistances <- function(sample, minBag, label) {

  sapply(minBag, function(y) {
    ifelse(y == sample,
      Inf, #In order not to choose its own
      sqrt( #Squared root because euclidean distance
        sum( #Summing distances between numeric and non numeric attributes
          sum( #For numeric attributes: square of the difference
            (D$dataset[sample,D$attributesIndexes[D$attributes[1:D$measures$num.inputs]=="numeric"]] - D$dataset[y,D$attributesIndexes[D$attributes[1:D$measures$num.inputs]=="numeric"]])^2
          ),
          sum( #For non numeric attributes: Value Difference Measure (VDM)
            sapply(D$attributesIndexes[D$attributes[1:D$measures$num.inputs]!="numeric"], function(x) {
              sum(
                abs((table((D$dataset[D$dataset[x] == D$dataset[sample,x],])[label])/(table(D$dataset[x])[[D$dataset[sample,x]]])) - (table((D$dataset[D$dataset[x] == D$dataset[y,x],])[label])/(table(D$dataset[x])[[D$dataset[y,x]]])))
              )
            })
          )
        )
      )
    )
  })

}



#' Auxiliary function used by MLSMOTE. Creates a synthetic sample based on values of attributes and labels of its neighbors
#'
#' @param sample Sample we are using as "template"
#' @param refNeigh Reference neighbor
#' @param neighbors Neighbors to take into account
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#'
#' @return A synthetic sample derived from the one passed as a parameter and its neighbors
#' @export
#' newSample(25,28,c(15,28,62),bibtex)
newSample <- function(sample, refNeigh, neighbors, D) {

  c(
    lapply(D$attributesIndexes, function(i) { #Attributes
      ifelse(D$attributes[[i]] %in% c("numeric", "Date"),
             D$dataset[sample,i] + D$dataset[refNeigh,i] - D$dataset[sample,i]*runif(1, 0, 1), #Numeric attributes
             tail(names(sort(table(D$dataset[neighbors, i]))), 1)) #Non numeric attributes
    }),
    sapply(lapply(D$dataset[c(sample, neighbors),D$labels$index], sum), function(x) { #Labels
      ifelse(x > ((length(neighbors)+1)/2), 1, 0)
    }),
    rep(NA, length(D$dataset) - D$measures$num.attributes) #Other measures like labelcount, SCUMBLE
  )

}



#' @title Synthetic oversampling of multilabel instances
#'
#' @description This function implements the MLSMOTE algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
#' whose aim is to identify instances with minoritary labels, and generate synthetic instances based on their neighbor instances.
#'
#' @source Francisco Charte, Antonio J. Rivera, María J. del Jesus, and Francisco Herrera. MLSMOTE: Approaching imbalanced multilabel learning through synthetic instance generation. Knowledge-Based Systems, 39:385–397, 2015. ISSN 0950-7051. doi:https://doi.org/10.1016/j.knosys.2015.07.019
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param k Number of neighbors to be considered when creatinga synthetic instance
#'
#' @return An mldr object containing the preprocessed multilabel dataset
#' @examples
#' library(mldr)
#' MLSMOTE(bibtex, 3)
#' @export
MLSMOTE <- function(D, k) {

  newSamples <- unlist(lapply(D$labels[D$labels$IRLbl > D$measures$meanIR,]$index, function(x) {
                                                                                    minBag <- as.numeric(rownames(D$dataset[D$dataset[x]==1,]))
                                                                                    lapply(minBag, function(y) {
                                                                                                    neighbors <- minBag[order(calculateDistances(y, minBag, x))[1:k+1]]
                                                                                                    refNeigh <- sample(neighbors, size=1)
                                                                                                    newSample(y, refNeigh, neighbors, D)
                                                                                                  })
                                                                                  }), recursive = FALSE)

  mldr_from_dataframe(rbind(D$dataset, lapply(setNames(as.data.frame(do.call(rbind, newSamples[-1])), names(D$dataset)), unlist)), D$labels$index, D$attributes, D$name)

}
