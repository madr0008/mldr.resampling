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

  minoritary <- unlist(.mldrApplyFun1(D$labels$freq, function(x) { ifelse(x<0.5,1,0) }, mc.cores=.numCores))
  numeric <- names(D$attributes[D$attributes == "numeric" | D$attributes == "numeric" ])
  factors <- names(D$attributes[D$attributes[D$attributesIndexes] != "numeric"])

  print("Part 1/3: Generating new instances")

  S <- .mldrApplyFun2(D$labels$index, function(l) {

    min <- which(D$dataset[[names(D$attributes)[l]]]==minoritary[l-D$measures$num.inputs])

    N <- D$measures$num.instances - 2*length(min)

    neighbors <- .mldrApplyFun1(min, function(i) {
      min[getNN(i, min, l, D, k)]
    }, mc.cores=.numCores)

    rNeighbors <- .mldrApplyFun1(min, function(i) {
      min[ceiling(which(unlist(neighbors)==i)/k)]
    }, mc.cores=.numCores)

    augMin <- min[lengths(rNeighbors) > 0]

    rNeighbors <- rNeighbors[.mldrApplyFun1(rNeighbors,length, mc.cores=.numCores)>0]

    q <- sample.int(length(augMin),N,replace=TRUE)

    aux <- as.data.frame(do.call(rbind, .mldrApplyFun1(q, function(i) {

      s <- augMin[i]
      r <- sample(rNeighbors[[i]],1)

      stats::setNames(c(unlist(.mldrApplyFun1(D$attributesIndexes, function(j) {

          if (D$attributes[[j]] %in% c("numeric", "Date")) {
            D$dataset[s,j] + (D$dataset[r,j] - D$dataset[s,j])*stats::runif(1, 0, 1) #Numeric attributes
          } else {
            sample(c(D$dataset[s,j], D$dataset[r,j]), size = 1) #Non numeric attributes
          }

      }, mc.cores=.numCores)), minoritary[l-D$measures$num.inputs]), names(D$attributes[c(D$attributesIndexes,l)]))

    }, mc.cores=.numCores)))

    aux[numeric] <- as.numeric(unlist(aux[numeric]))
    aux[factors] <- .mldrApplyFun1(factors, function(i) { factor(aux[[i]], levels=unique(D$dataset[[i]])) }, mc.cores=.numCores)
    aux

  }, mc.cores=.numCores)

  print("Part 2/3: Training SVMs")

  #Train SVM
  M <- .mldrApplyFun2(D$labels$index, function(l) {

    d <- D$dataset[,c(D$attributesIndexes,l)]
    d[factors] <- .mldrApplyFun1(factors, function(i) { factor(d[[i]], levels=unique(D$dataset[[i]])) }, mc.cores=.numCores)

    O <- rbind(D$dataset[,c(D$attributesIndexes,l)], S[[l-D$measures$num.inputs]])
    names(O) <- c(names(D$attributes[D$attributesIndexes]),"Class")
    O$Class <- as.factor(O$Class)

    e1071::svm(Class ~ ., data = O)

  }, mc.cores=.numCores)

  print("Part 3/3: Predicting labels")

  toRet <- as.data.frame(data.table::rbindlist(list(D$dataset, data.table::rbindlist(.mldrApplyFun1(D$labels$index, function(j) {

    S[[j - D$measures$num.inputs]][[names(D$dataset)[j]]] <- as.numeric(S[[j - D$measures$num.inputs]][[names(D$dataset)[j]]])

    cbind(S[[j - D$measures$num.inputs]], stats::setNames(.mldrApplyFun1(D$labels$index[!D$labels$index %in% j], function(l) {

      as.numeric(stats::predict(M[[l - D$measures$num.inputs]],(S[[j - D$measures$num.inputs]])[D$attributesIndexes]))

    }, mc.cores=.numCores), names(D$dataset)[D$labels$index[!D$labels$index %in% j]]))


  }, mc.cores=.numCores), fill=TRUE, use.names=TRUE)), fill=TRUE, use.names=TRUE))

  toRet[] <- .mldrApplyFun1(toRet, function(x) { if (is.factor(x)) as.character(x) else {x} }, mc.cores=.numCores)

  mldr::mldr_from_dataframe(toRet, D$labels$index, D$attributes, D$name)

}
