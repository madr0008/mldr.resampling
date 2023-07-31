#' @title Reverse-nearest neighborhood based oversampling for imbalanced, multi-label datasets
#'
#' @description This function implements an algorithm that uses the concept of reverse nearest neighbors, in order to create new instances for each label. Then, several radial SVMs, one for each label, are trained in order to predict each label of the synthetic instances.
#'
#' @source Sadhukhan, P., & Palit, S. (2019). Reverse-nearest neighborhood based oversampling for imbalanced, multi-label datasets. Pattern Recognition Letters, 125, 813-820
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param k Number of neighbors to be considered when creating a synthetic instance
#' @param tableVDM Dataframe object containing previous calculations for faster processing. If it is empty, the algorithm will be slower
#'
#' @return A mld object containing the preprocessed multilabel dataset
#'
#' @export
MLRkNNOS <- function(D, k, tableVDM=NULL) {

  minoritary <- unlist(mldr.resampling.env$.mldrApplyFun1(D$labels$freq, function(x) { ifelse(x<0.5,1,0) }, mc.cores=mldr.resampling.env$.numCores))
  numeric <- names(D$attributes[D$attributes == "numeric" | D$attributes == "numeric" ])
  factors <- names(D$attributes[D$attributes[D$attributesIndexes] != "numeric"])

  message("Part 1/3: Generating new instances")

  S <- mldr.resampling.env$.mldrApplyFun2(D$labels$index, function(l) {

    min <- which(D$dataset[[names(D$attributes)[l]]]==minoritary[l-D$measures$num.inputs])

    N <- D$measures$num.instances - 2*length(min)

    neighbors <- mldr.resampling.env$.mldrApplyFun1(min, function(i) {
      min[getNN(i, min, l, D, tableVDM)[1:k+1]]
    }, mc.cores=mldr.resampling.env$.numCores)

    rNeighbors <- mldr.resampling.env$.mldrApplyFun1(min, function(i) {
      min[ceiling(which(unlist(neighbors)==i)/k)]
    }, mc.cores=mldr.resampling.env$.numCores)

    augMin <- min[lengths(rNeighbors) > 0]

    rNeighbors <- rNeighbors[mldr.resampling.env$.mldrApplyFun1(rNeighbors,length, mc.cores=mldr.resampling.env$.numCores)>0]

    q <- sample.int(length(augMin),N,replace=TRUE)

    aux <- as.data.frame(do.call(rbind, mldr.resampling.env$.mldrApplyFun1(q, function(i) {

      s <- augMin[i]
      r <- sample(rNeighbors[[i]],1)

      stats::setNames(c(unlist(mldr.resampling.env$.mldrApplyFun1(D$attributesIndexes, function(j) {

          if (D$attributes[[j]] %in% c("numeric", "Date")) {
            D$dataset[s,j] + (D$dataset[r,j] - D$dataset[s,j])*stats::runif(1, 0, 1) #Numeric attributes
          } else {
            sample(c(D$dataset[s,j], D$dataset[r,j]), size = 1) #Non numeric attributes
          }

      }, mc.cores=mldr.resampling.env$.numCores)), minoritary[l-D$measures$num.inputs]), names(D$attributes[c(D$attributesIndexes,l)]))

    }, mc.cores=mldr.resampling.env$.numCores)))

    aux[numeric] <- as.numeric(unlist(aux[numeric]))
    aux[factors] <- mldr.resampling.env$.mldrApplyFun1(factors, function(i) { factor(aux[[i]], levels=unique(D$dataset[[i]])) }, mc.cores=mldr.resampling.env$.numCores)
    aux

  }, mc.cores=mldr.resampling.env$.numCores, parL=c("D", "k", "tableVDM", "minoritary", "numeric", "factors"), parEnv=environment())

  message("Part 2/3: Training SVMs")

  #Train SVM
  M <- mldr.resampling.env$.mldrApplyFun2(D$labels$index, function(l) {

    d <- D$dataset[,c(D$attributesIndexes,l)]
    d[factors] <- mldr.resampling.env$.mldrApplyFun1(factors, function(i) { factor(d[[i]], levels=unique(D$dataset[[i]])) }, mc.cores=mldr.resampling.env$.numCores)

    O <- rbind(D$dataset[,c(D$attributesIndexes,l)], S[[l-D$measures$num.inputs]])
    names(O) <- c(names(D$attributes[D$attributesIndexes]),"Class")
    O$Class <- as.factor(O$Class)

    e1071::svm(Class ~ ., data = O)

  }, mc.cores=mldr.resampling.env$.numCores, parL=c("D", "factors", "S"), parEnv=environment())

  message("Part 3/3: Predicting labels")

  toRet <- as.data.frame(data.table::rbindlist(list(D$dataset, data.table::rbindlist(mldr.resampling.env$.mldrApplyFun1(D$labels$index, function(j) {

    S[[j - D$measures$num.inputs]][[names(D$dataset)[j]]] <- as.numeric(S[[j - D$measures$num.inputs]][[names(D$dataset)[j]]])

    cbind(S[[j - D$measures$num.inputs]], stats::setNames(mldr.resampling.env$.mldrApplyFun1(D$labels$index[!D$labels$index %in% j], function(l) {

      as.numeric(stats::predict(M[[l - D$measures$num.inputs]],(S[[j - D$measures$num.inputs]])[D$attributesIndexes]))

    }, mc.cores=mldr.resampling.env$.numCores), names(D$dataset)[D$labels$index[!D$labels$index %in% j]]))


  }, mc.cores=mldr.resampling.env$.numCores), fill=TRUE, use.names=TRUE)), fill=TRUE, use.names=TRUE))

  toRet[] <- mldr.resampling.env$.mldrApplyFun1(toRet, function(x) { if (is.factor(x)) as.character(x) else {x} }, mc.cores=mldr.resampling.env$.numCores)

  mldr::mldr_from_dataframe(toRet, D$labels$index, D$attributes, D$name)

}
