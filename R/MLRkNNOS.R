#' @title Reverse-nearest neighborhood based oversampling for imbalanced, multi-label datasets
#'
#' @description This function implements an algorithm that uses the concept of reverse nearest neighbors, in order to create new instances for each label. Then, several radial SVMs, one for each label, are trained in order to predict each label of the synthetic instances.
#'
#' @source Sadhukhan, P., & Palit, S. (2019). Reverse-nearest neighborhood based oversampling for imbalanced, multi-label datasets. Pattern Recognition Letters, 125, 813-820
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param k Number of neighbors to be considered when creating a synthetic instance
#'
#' @return A mld object containing the preprocessed multilabel dataset
#' @examples
#' \dontrun{
#' library(mldr)
#' MLRkNNOS(bibtex, 3)
#' }
#' @export
MLRkNNOS <- function(D, k) {

  minoritary <- unlist(lapply(D$labels$freq, function(x) ifelse(x<0.5,1,0)))
  numeric <- names(D$attributes[D$attributes == "numeric" | D$attributes == "numeric" ])
  factors <- names(D$attributes[D$attributes[D$attributesIndexes] != "numeric"])

  print("Part 1/3: Generating new instances")

  S <- pbapply::pblapply(D$labels$index, function(l) {

    min <- which(D$dataset[[names(D$attributes)[l]]]==minoritary[l-D$measures$num.inputs])

    N <- D$measures$num.instances - 2*length(min)

    neighbors <- lapply(min, function(i) {
      min[getNN(i, min, l, D, k)]
    })

    rNeighbors <- lapply(min, function(i) {
      min[ceiling(which(unlist(neighbors)==i)/k)]
    })

    augMin <- min[lengths(rNeighbors) > 0]

    rNeighbors <- rNeighbors[lapply(rNeighbors,length)>0]

    q <- sample.int(length(augMin),N,replace=TRUE)

    aux <- as.data.frame(do.call(rbind, lapply(q, function(i) {

      s <- augMin[i]
      r <- sample(rNeighbors[[i]],1)

      stats::setNames(c(unlist(lapply(D$attributesIndexes, function(j) {

          if (D$attributes[[j]] %in% c("numeric", "Date")) {
            D$dataset[s,j] + (D$dataset[r,j] - D$dataset[s,j])*stats::runif(1, 0, 1) #Numeric attributes. Falta multiplicar por dist euclidea
          } else {
            sample(c(D$dataset[s,j], D$dataset[r,j]), size = 1) #Non numeric attributes
          }

      })), minoritary[l-D$measures$num.inputs]), names(D$attributes[c(D$attributesIndexes,l)]))

    })))

    aux[numeric] <- as.numeric(unlist(aux[numeric]))
    aux[factors] <- lapply(factors, function(i) { factor(aux[[i]], levels=unique(D$dataset[[i]])) })
    aux

  })

  print("Part 2/3: Training SVMs")

  #Train SVM
  M <- pbapply::pblapply(D$labels$index, function(l) {

    d <- D$dataset[,c(D$attributesIndexes,l)]
    d[factors] <- lapply(factors, function(i) { factor(d[[i]], levels=unique(D$dataset[[i]])) })

    O <- rbind(D$dataset[,c(D$attributesIndexes,l)], S[[l-D$measures$num.inputs]])
    names(O) <- c(names(D$attributes[D$attributesIndexes]),"Class")
    O$Class <- as.factor(O$Class)

    e1071::svm(Class ~ ., data = O)

  })

  print("Part 3/3: Predicting labels")

  toRet <- as.data.frame(data.table::rbindlist(list(D$dataset, data.table::rbindlist(lapply(D$labels$index, function(j) {

    S[[j - D$measures$num.inputs]][[names(D$dataset)[j]]] <- as.numeric(S[[j - D$measures$num.inputs]][[names(D$dataset)[j]]])

    cbind(S[[j - D$measures$num.inputs]], stats::setNames(lapply(D$labels$index[!D$labels$index %in% j], function(l) {

      as.numeric(stats::predict(M[[l - D$measures$num.inputs]],(S[[j - D$measures$num.inputs]])[D$attributesIndexes]))

    }), names(D$dataset)[D$labels$index[!D$labels$index %in% j]]))


  }), fill=TRUE, use.names=TRUE)), fill=TRUE, use.names=TRUE))

  toRet[] <- lapply(toRet, function(x) if (is.factor(x)) as.character(x) else {x})

  mldr::mldr_from_dataframe(toRet, D$labels$index, D$attributes, D$name)

}
