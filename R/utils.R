#' Auxiliary function used to calculate the distances between an instance and the ones with a specific active label. Euclidean distance is calculated for numeric attributes, and VDM for non numeric ones.
#'
#' @param sample Index of the sample whose distances to other samples we want to know
#' @param rest Indexes of the samples to which we will calculate the distance
#' @param label Label that must be active
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param tableVDM Dataframe object containing previous calculations for faster processing. If it is empty, the algorithm will be slower
#'
#' @return A list with the distance to the rest of samples
calculateDistances <- function(sample, rest, label, D, tableVDM=NULL) {

  unlist(mldr.resampling.env$.mldrApplyFun1(rest, function(y) {
    ifelse(y == sample,
           Inf, #In order not to choose its own
           sqrt( #Square root because euclidean distance
             sum( #Summing distances between numeric and non numeric attributes
               sum((D$dataset[sample,D$attributesIndexes[D$attributes[1:D$measures$num.inputs]=="numeric"]] - D$dataset[y,D$attributesIndexes[D$attributes[1:D$measures$num.inputs]=="numeric"]])^2), #For numeric attributes: square of the difference
               ifelse(sum(D$attributes[1:D$measures$num.inputs]!="numeric") > 0, vdm(D,sample,y,label,tableVDM), 0) #For non numeric attributes: Value Difference Metric (VDM)
             )
           )
    )
  }, mc.cores=mldr.resampling.env$.numCores))

}



#' Auxiliary function used to calculate the Value Difference Metric (VDM) between two instances considering their non numeric attributes
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param sample Index of the first sample
#' @param y Index of the second sample
#' @param label Label that will be considered in calculations
#' @param tableVDM Dataframe object containing previous calculations for faster processing. If it is empty, the algorithm will be slower
#'
#' @return A value for the distance
vdm <- function(D, sample, y, label, tableVDM=NULL) {

  if (is.null(tableVDM)) {

    sum(
      unlist(mldr.resampling.env$.mldrApplyFun1(D$attributesIndexes[1:D$measures$num.inputs][D$attributes[1:D$measures$num.inputs]!="numeric"], function(x) {
        table1 <- table((D$dataset[D$dataset[x] == D$dataset[sample,x],])[label])/(table(D$dataset[x])[[D$dataset[sample,x]]])
        table2 <- table((D$dataset[D$dataset[x] == D$dataset[y,x],])[label])/(table(D$dataset[x])[[D$dataset[y,x]]])
        sum(
          abs(ifelse(length(table1 == 1), ifelse(names(table1) == "0", stats::setNames(c(table1, 0), c("0","1")), stats::setNames(c(0, table1), c("0","1"))), table1) - ifelse(length(table2 == 1), ifelse(names(table2) == "0", stats::setNames(c(table2, 0), c("0","1")), stats::setNames(c(0, table2), c("0","1"))), table2))
        )
      }, mc.cores=mldr.resampling.env$.numCores))
    )

  } else {

    sum( #For non numeric attributes: Value Difference Metric (VDM)
      unlist(mldr.resampling.env$.mldrApplyFun1(tableVDM$attribute, function(x) {
        table1 <- (tableVDM[tableVDM$attribute==x,]$tablePerLabel[[1]][[match(D$dataset[sample,x],names(unlist(tableVDM[tableVDM$attribute==x,]$generalTable)))]][[label-D$measures$num.inputs]])/(unlist(tableVDM[tableVDM$attribute==x,]$generalTable)[[D$dataset[sample,x]]])
        table2 <- (tableVDM[tableVDM$attribute==x,]$tablePerLabel[[1]][[match(D$dataset[y,x],names(unlist(tableVDM[tableVDM$attribute==x,]$generalTable)))]][[label-D$measures$num.inputs]])/(unlist(tableVDM[tableVDM$attribute==x,]$generalTable)[[D$dataset[y,x]]])
        sum(
          abs(ifelse(length(table1 == 1), ifelse(names(table1) == "0", stats::setNames(c(table1, 0), c("0","1")), stats::setNames(c(0, table1), c("0","1"))), table1) - ifelse(length(table2 == 1), ifelse(names(table2) == "0", stats::setNames(c(table2, 0), c("0","1")), stats::setNames(c(0, table2), c("0","1"))), table2))
        )
      }, mc.cores=mldr.resampling.env$.numCores))
    )

  }



}



#' Auxiliary function used to calculate an auxiliary table to make VDM calculation faster
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#'
#' @return A dataframe with tables, useful for VDM calculation
calculateTableVDM <- function(D) {

  attribute <- D$attributesIndexes[1:D$measures$num.inputs][D$attributes[1:D$measures$num.inputs]!="numeric"]
  df <- data.frame(attribute)

  df$generalTable <- mldr.resampling.env$.mldrApplyFun1(attribute, function(x) { table(D$dataset[x]) },mc.cores=mldr.resampling.env$.numCores)

  df$tablePerLabel <- mldr.resampling.env$.mldrApplyFun2(attribute, function(x) {

  mldr.resampling.env$.mldrApplyFun1(unlist(unname(unique(D$dataset[x]))), function(y) {

    mldr.resampling.env$.mldrApplyFun1(D$labels$index, function(l) {

      table(D$dataset[D$dataset[x]==y,l])

      }, mc.cores=mldr.resampling.env$.numCores)

    }, mc.cores=mldr.resampling.env$.numCores)

  }, mc.cores=mldr.resampling.env$.numCores)

  df

}



#' Auxiliary function used to compute the neighbors of an instance
#'
#' @param sample Index of the sample whose neighbors we want to know
#' @param rest Indexes of the samples among which we will search
#' @param label Label that must be active, in order to calculate the distances
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param tableVDM Dataframe object containing previous calculations for faster processing. If it is empty, the algorithm will be slower
#'
#' @return A vector with the indexes inside rest of the neighbors
getNN <- function(sample, rest, label, D, tableVDM=NULL) {

  order(calculateDistances(sample, rest, label, D, tableVDM))

}



#' Auxiliary function used by MLSMOTE. Creates a synthetic sample based on values of attributes and labels of its neighbors
#'
#' @param seedInstance Sample we are using as "template"
#' @param refNeigh Reference neighbor
#' @param neighbors Neighbors to take into account
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#'
#' @return A synthetic sample derived from the one passed as a parameter and its neighbors
newSample <- function(seedInstance, refNeigh, neighbors, D) {

  c(
    mldr.resampling.env$.mldrApplyFun1(D$attributesIndexes[1:D$measures$num.inputs], function(i) { #Attributes
      ifelse(D$attributes[[i]] %in% c("numeric", "Date"),
             D$dataset[seedInstance,i] + (D$dataset[refNeigh,i] - D$dataset[seedInstance,i])*stats::runif(1, 0, 1), #Numeric attributes
             utils::tail(names(sort(table(D$dataset[neighbors, i]))), 1)) #Non numeric attributes
    }, mc.cores=mldr.resampling.env$.numCores),
    unlist(mldr.resampling.env$.mldrApplyFun1(mldr.resampling.env$.mldrApplyFun1(D$dataset[c(seedInstance, neighbors),D$labels$index], sum, mc.cores=mldr.resampling.env$.numCores), function(x) { #Labels
      ifelse(x > ((length(neighbors)+1)/2), 1, 0)
    }, mc.cores=mldr.resampling.env$.numCores),
  ))

}



#' Auxiliary function used by MLeNN. Computes the Hamming Distance between two instances
#'
#' @param x Index of sample 1
#' @param y Index of sample 2
#' @param D mld \code{mldr} object in which the instances are located
#'
#' @return The Hamming Distance between the instances
adjustedHammingDist <- function(x,y,D) {
  length(which(D$dataset[x,D$labels$index] != D$dataset[y,D$labels$index])) / (sum(D$dataset[x,D$labels$index]) + sum(D$dataset[y,D$labels$index]))
}



#' Auxiliary function used by MLSOL. Categorizes each pair instance-label of the dataset with a type
#'
#' @param C List of vectors with one value for each pair instance-label
#' @param neighbors Structure with the k nearest neighbors of each instance of the dataset
#' @param k Number of neighbors to be considered for each instance
#' @param minoritary Vector with the minoritary value of each label (normally, 1)
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#'
#' @return A synthetic sample derived from the one passed as a parameter and its neighbors
initTypes <- function(C, neighbors, k, minoritary, D, d) {

  t <- mldr.resampling.env$.mldrApplyFun2(c(1:length(d)), function(i) {
    unlist(mldr.resampling.env$.mldrApplyFun1(D$labels$index, function(j) {
      if (D$dataset[d[i],j] != minoritary[j - D$measures$num.inputs]) {
        5
      } else {
        if (C[[i]][j - D$measures$num.inputs] < 1/3) {
          1
        } else if (C[[i]][j - D$measures$num.inputs] < 2/3) {
          2
        } else if (C[[i]][j - D$measures$num.inputs] < 1) {
          3
        } else {
          4
        }
      }
    }, mc.cores=mldr.resampling.env$.numCores))
  }, mc.cores=mldr.resampling.env$.numCores)

  change <- TRUE
  while (change) {
    change <- FALSE
    for (i in c(1:length(d))) {
      for (j in 1:D$measures$num.labels) {
        if (t[[i]][j] == 3) {
          for (m in neighbors[[i]]) {
            if (t[[m]][j] == 1 | t[[m]][j] == 2) {
              t[[i]][j] <- 2
              change <- TRUE
              break
            }
          }
        }
      }
    }
  }

  t

}



#' Auxiliary function used by MLSOL. Creates a synthetic sample based on two other samples, taking into account their types
#'
#' @param seedInstance Index of the sample we are using as "template"
#' @param refNeigh Index of the reference neighbor
#' @param t types of the instances
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#'
#' @return A synthetic sample derived from the one passed as a parameter and its neighbors
generateInstanceMLSOL <- function (seedInstance, refNeigh, t, D) {

  s <- seedInstance
  r <- refNeigh

  attributes <- mldr.resampling.env$.mldrApplyFun1(D$attributesIndexes[1:D$measures$num.inputs], function(i) { #Attributes
    ifelse(D$attributes[[i]] %in% c("numeric", "Date"),
           D$dataset[s,i] + (D$dataset[r,i] - D$dataset[s,i])*stats::runif(1, 0, 1), #Numeric attributes
           sample(c(D$dataset[s,i], D$dataset[r,i]), size = 1)) #Non numeric attributes
  }, mc.cores=mldr.resampling.env$.numCores)

  #Calculate distances between attributes
  numericAttributes <- D$attributesIndexes[1:D$measures$num.inputs][D$attributes[1:D$measures$num.inputs] %in% c("numeric", "Date")]
  categoricAttributes <- D$attributesIndexes[1:D$measures$num.inputs][!(D$attributes[1:D$measures$num.inputs] %in% c("numeric", "Date"))]

  d_s <- sum(
           ifelse(length(numericAttributes) > 0, sum(attributes[numericAttributes] - D$dataset[s,numericAttributes])^2, 0),
           sum(unlist(mldr.resampling.env$.mldrApplyFun1(categoricAttributes, function(j) {ifelse(attributes[[j]] == D$dataset[s,j], 0, 1) })))
         )

  d_r <- sum(
    ifelse(length(numericAttributes) > 0, sum(attributes[numericAttributes] - D$dataset[r,numericAttributes])^2, 0),
    sum(unlist(mldr.resampling.env$.mldrApplyFun1(categoricAttributes, function(j) {ifelse(attributes[[j]] == D$dataset[r,j], 0, 1) })))
  )

  cd <- d_s / (d_s + d_r)

  labels <- c()

  for (j in D$labels$index) {
    if (D$dataset[s,j] == D$dataset[r,j]) {
      labels <- c(labels, D$dataset[s,j])
    } else {
      if (t[[s]][j - D$measures$num.inputs] == 5) {
        aux <- s
        s <- r
        r <- aux
        cd <- 1 - cd
      }
      theta <- 0.5
      if (t[[s]][j - D$measures$num.inputs] == 2) {
        theta <- 0.75
      } else if (t[[s]][j - D$measures$num.inputs] == 3) {
        theta <- 1.00001
      } else if (t[[s]][j - D$measures$num.inputs] == 4) {
        theta <- -0.00001
      }
      if (cd <= theta) {
        labels <- c(labels, D$dataset[s,j])
      } else {
        labels <- c(labels, D$dataset[r,j])
      }
    }
  }

  c(attributes, labels)

}



#' Auxiliary function used by MLSOL and MLUL. Computes the kNN of every instance in a dataset
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param tableVDM Dataframe object containing previous calculations for faster processing. If it is empty, the algorithm will be slower
#'
#' @return A list of vectors with the indexes of the neighbors for each instance
getAllNeighbors <- function(D, d, tableVDM=NULL) {

  mldr.resampling.env$.mldrApplyFun2(d, function(i) {
    activeLabels <- D$labels[which(D$dataset[i,D$labels$index] %in% 1),1]
    if (length(activeLabels) > 0) {
      getNN(i, d, ifelse(length(activeLabels)==1,activeLabels,sample(activeLabels,1)), D, tableVDM)
    }
  }, mc.cores=mldr.resampling.env$.numCores)

}



#' Auxiliary function used by MLeNN and MLTL. Gets the kNN of every instance in a dataset, when compared to some of the rest
#'
#' @param neighbors Structure with all the neighbors in the dataset, regardless of which ones to be compared
#' @param d Vector with the instances of the dataset which are going to be compared
#' @param k Number of neighbors to be retrieved
#'
#' @return A list of vectors with the indexes of the neighbors for each instance
getAllNeighbors2 <- function(neighbors, d, k) {

  mldr.resampling.env$.mldrApplyFun2(neighbors, function(x) {
    intersect(x, d)[1:k+1]
  }, mc.cores=mldr.resampling.env$.numCores)

}



#' Auxiliary function used by MLSOL and MLUL. For each instance in the dataset, we compute, for each label, the proportion of neighbors having an opposite class with respect to the proper instance
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param neighbors Structure with the neighbors of every instance in the dataset
#' @param k Number of neighbors taken into account for each instance
#'
#' @return A structure with the proportion of neighbors having an opposite class with respect to an instance and label
getC <- function(D, d, neighbors, k) {

  mldr.resampling.env$.mldrApplyFun2(d, function(i) {
    unlist(mldr.resampling.env$.mldrApplyFun1(D$labels$index, function(j) {
      (1/k) * sum(unlist(mldr.resampling.env$.mldrApplyFun1(neighbors[i], function(m) {
        ifelse(D$dataset[i,j]==D$dataset[m,j],0,1)
      }, mc.cores=mldr.resampling.env$.numCores)))
    }, mc.cores=mldr.resampling.env$.numCores))
  }, mc.cores=mldr.resampling.env$.numCores)

}



#' Auxiliary function used by MLSOL and MLUL. For non outlier instances, it aggregates the values of C, taking into account the global class imbalance
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param C Structure with the proportion of neighbors having an opposite class with respect to an instance and label
#' @param minoritary Vector with the minoritary class of each label (normally, 1)
#'
#' @return A structure with the proportion of neighbors having an opposite class with respect to an instance and label, normalized by the global class imbalance
getS <- function(D, d, C, minoritary) {

  mldr.resampling.env$.mldrApplyFun2(c(1:length(d)), function(i) {
    unlist(mldr.resampling.env$.mldrApplyFun1(D$labels$index, function(j) {
      if ((D$dataset[d[i],j]==minoritary[j - D$measures$num.inputs]) & (C[[i]][j - D$measures$num.inputs]<1)) {
        numerator <- C[[i]][j - D$measures$num.inputs]
        denominator <- sum(unlist(mldr.resampling.env$.mldrApplyFun1(c(1:length(d)), function(x) {
          ifelse(D$dataset[d[x],j]==minoritary[j - D$measures$num.inputs],C[[x]][j - D$measures$num.inputs],0)*ifelse(C[[x]][j - D$measures$num.inputs]<1,1,0)
        }, mc.cores=mldr.resampling.env$.numCores)))
        ifelse(denominator!=0,numerator/denominator,-1)
      } else {
        -1
      }
    }, mc.cores=mldr.resampling.env$.numCores))
  }, mc.cores=mldr.resampling.env$.numCores)

}



#' Auxiliary function used by MLSOL and MLUL. For non outlier instances, it aggregates the values of S for each label
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param S Structure with the proportion of neighbors having an opposite class with respect to an instance and label, normalized by the global class imbalance
#'
#' @return A vector of weights to be considered when oversampling for each instance
getW <- function(D, S) {

  mldr.resampling.env$.mldrApplyFun2(c(1:length(S)), function(i) {
    sum(S[[i]][!S[[i]] %in% -1])
  }, mc.cores=mldr.resampling.env$.numCores)

}



#' Auxiliary function used by MLUL. For each instance in the dataset, given the neighbors structure, we compute its reverse nearest neighbors
#'
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param neighbors Structure with the neighbors of every instance in the dataset
#' @param k Number of neighbors to be considered
#'
#' @return A list of vectors with the indexes of the reverse nearest neighbors of every instance in the dataset
getAllReverseNeighbors <- function(d, neighbors, k) {

  mldr.resampling.env$.mldrApplyFun2(d, function(i) {
    ceiling(which(unlist(neighbors)==i)/k)
  }, mc.cores=mldr.resampling.env$.numCores)

}



#' Auxiliary function used by MLUL. It computes the influence of each instance with respect to its reverse neighbors
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param rNeighbors Structure with the reverse nearest neighbors of each instance of the dataset
#' @param S Structure with the proportion of neighbors having an opposite class with respect to an instance and label, normalized by the global class imbalance
#'
#' @return A list of values of influence for each instance with respect to its reverse neighbors
getU <- function(D, d, rNeighbors, S) {

  mldr.resampling.env$.mldrApplyFun2(c(1:length(d)), function(i) {

    sum(unlist(mldr.resampling.env$.mldrApplyFun1(D$labels$index, function(j) {

        ifelse(length(rNeighbors[[i]]) > 0,

        sum(unlist(mldr.resampling.env$.mldrApplyFun1(rNeighbors[[i]], function(m) {

          ifelse(S[[m]][j - D$measures$num.inputs]==-1,0,ifelse(D$dataset[i,j]==D$dataset[m,j],1,-1)*S[[m]][j - D$measures$num.inputs])

        }, mc.cores=mldr.resampling.env$.numCores))) / length(rNeighbors[i]), 0)

    }, mc.cores=mldr.resampling.env$.numCores)))

  }, mc.cores=mldr.resampling.env$.numCores)

}



#' Auxiliary function used by MLUL. It calculates, for each instance, how important it is in the dataset
#'
#' @param w List of weights for each instance
#' @param u List of influences in reverse neighbors for each instance
#'
#' @return A list with the values of importance of each instance in the dataset
getV <- function(w, u) {

  v <- mldr.resampling.env$.mldrApplyFun2(c(1:length(w)), function(i) {
         w[[i]] + u[[i]]
       }, mc.cores=mldr.resampling.env$.numCores)

  minimum <- min(unlist(v))

  mldr.resampling.env$.mldrApplyFun1(v, function(x) {
    x - minimum
  }, mc.cores=mldr.resampling.env$.numCores)

}



#' Auxiliary function used by resample. It executes an algorithm, given as a string, and stores the resulting MLD in a arff file
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param a String with the name of the algorithm to be applied.
#' @param P Percentage in which the original dataset is increased/decreased (if required by the algorithm)
#' @param k Number of neighbors taken into account for each instance (if required by the algorithm)
#' @param TH Threshold for the Hamming Distance in order to consider an instance different to another one (if required by the algorithm)
#' @param outputDirectory Route with the directory where the generated ARFF file will be stored
#' @param neighbors Structure with all instances and neighbors in the dataset, useful in MLSOL and MLUL
#' @param neighbors2 Structure with some instances and neighbors in the dataset, useful in MLeNN and MLTL
#' @param tableVDM Dataframe object containing previous calculations for faster processing. If it is empty, the algorithm will be slower
#'
#' @return Time (in seconds) taken to execute the algorithm (NULL if no algorithm was executed)
executeAlgorithm <- function(D, a, P, k, TH, outputDirectory, neighbors, neighbors2, tableVDM) {

  if (!(a %in% c("LPROS", "LPRUS", "MLROS", "MLRUS", "MLRkNNOS", "MLSMOTE", "MLSOL", "MLUL", "MLeNN", "MLTL", "REMEDIAL"))) {
    message(paste("Error: There is no algorithm named", a))
    NULL
  } else {

    f <- get(a)
    if (a %in% c("LPROS", "LPRUS", "MLROS", "MLRUS")) {
      name <- paste(D$name, a, "P", P, sep = "_")
      message(paste("Running",a,"on",D$name,"with P =",P))
      startTime <- Sys.time()
      d <- f(D, P)
      endTime <- Sys.time()
    } else if (a %in% c("MLRkNNOS", "MLSMOTE")) {
      name <- paste(D$name, a, "k", k, sep = "_")
      message(paste("Running",a,"on",D$name,"with k =",k))
      startTime <- Sys.time()
      d <- f(D, k, tableVDM)
      endTime <- Sys.time()
    } else if (a %in% c("MLSOL", "MLUL")) {
      name <- paste(D$name, a, "P", P, "k", k, sep = "_")
      message(paste("Running",a,"on",D$name,"with P =",P,"and k =",k))
      startTime <- Sys.time()
      d <- f(D, P, k, neighbors, tableVDM)
      endTime <- Sys.time()
    } else if (a == "MLeNN") {
      name <- paste(D$name, a, "TH", TH, "k", k, sep = "_")
      message(paste("Running",a,"on",D$name,"with TH =",TH,"and k =",k))
      startTime <- Sys.time()
      d <- f(D, TH, k, neighbors2, tableVDM)
      endTime <- Sys.time()
    } else if (a == "MLTL") {
      name <- paste(D$name, a, "TH", TH, sep = "_")
      message(paste("Running",a,"on",D$name,"with TH =",TH))
      startTime <- Sys.time()
      d <- f(D, TH, neighbors2, tableVDM)
      endTime <- Sys.time()
    } else { #REMEDIAL
      name <- a
      message(paste("Running",a,"on",D$name))
      startTime <- Sys.time()
      d <- f(D)
      endTime <- Sys.time()
    }

    time <- as.numeric(endTime - startTime, units="secs")

    message(paste("Time taken (in seconds):",time))

    mldr::write_arff(d, paste(outputDirectory, name, sep="/"))

    time

  }

}



#' Interface function of the package. It executes one or several algorithms, given as strings, and stores the resulting MLDs in arff files
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param algorithms String, or string vector, with the name(s) of the algorithm(s) to be applied.
#' @param P Percentage in which the original dataset is increased/decreased, if required by the algorithm(s). Defaults to 25
#' @param k Number of neighbors taken into account for each instance, if required by the algorithm(s). Defaults to 3
#' @param TH Threshold for the Hamming Distance in order to consider an instance different to another one, if required by the algorithm(s). Defaults to 0.5
#' @param params Dataframe with 4 columns: name of the algorithm, P, k and TH, in that order, to execute several algorithms with different values for their parameters
#' @param outputDirectory Route with the directory where generated ARFF files will be stored. Defaults to a temporary directory
#'
#' @return Dataframe with times (in seconds) taken in to execute each algorithm
#'
#' @examples
#' library(mldr)
#' library(mldr.resampling)
#' resample(birds, "LPROS", P=25)
#' resample(birds, c("LPROS", "LPRUS"), P=30)
#' @export
resample <- function(D, algorithms, P=25, k=3, TH=0.5, params, outputDirectory=tempdir()) {

  times <- data.frame(matrix(nrow = 0, ncol = 2))
  colnames(times) <- c("algorithm", "time")

  if (missing(D)) {

    message("Please, provide a mld object as the original dataset")
    NULL

  } else {

    if (missing(algorithms)) {

      if (missing(params)) {

        message("Please, specify the dataset and algorithms to be applied, either with the parameter algorithms or with params")
        NULL

      } else {

        tableVDM <- NULL
        timeTableVDM <- 0
        neighbors <- NULL
        neighbors2 <- NULL
        timeNeighbors <- 0

        message(paste("Calculating structures for dataset", D$name, ", if necessary. Once this is done, algorithms will be applied faster"))

        if (length(intersect(algorithms, c("MLSMOTE","MLeNN","MLTL","MLSOL","MLUL"))) > 0) {

          message(paste("Calculating VDM table for dataset", D$name))
          startTime <- Sys.time()
          tableVDM <- calculateTableVDM(D)
          endTime <- Sys.time()
          timeTableVDM <- as.numeric(endTime - startTime, units="secs")
          message(paste("Time taken (in seconds):",timeTableVDM))

        }

        if (length(intersect(algorithms, c("MLeNN", "MLTL","MLSOL","MLUL"))) > 0) {

          message(paste("Calculating neighbors structure for dataset", D$name))
          startTime <- Sys.time()
          neighbors <- getAllNeighbors(D, c(1:D$measures$num.instances)[D$dataset$.labelcount > 0], tableVDM)
          neighbors2 <- getAllNeighbors2(neighbors, unique(unlist(mldr.resampling.env$.mldrApplyFun1(D$labels[D$labels$IRLbl < D$measures$meanIR,]$index, function(x) { c(1:D$measures$num.instances)[D$dataset[x]==1] }, mc.cores=mldr.resampling.env$.numCores))), max(params[params[,1] %in% c("MLeNN","MLTL"),3]))
          neighbors <- mldr.resampling.env$.mldrApplyFun1(neighbors, function(x) { x[1:(max(params[params[,1] %in% c("MLSOL","MLUL"),3])+1)] }, mc.cores=mldr.resampling.env$.numCores)
          endTime <- Sys.time()
          timeNeighbors <- as.numeric(endTime - startTime, units="secs")
          message(paste("Time taken (in seconds):",timeNeighbors))

        }

        for(i in 1:nrow(params)) {

          time <- executeAlgorithm(D,params[i,1],params[i,2],params[i,3],params[i,4], neighbors, neighbors2)

          if (!is.null(time)) {
            times[nrow(times) + 1,] <- c(a,time)
          }

        }

        #Add structure generation times
        for (i in rownames(times[times$algorithm %in% c("MLSOL","MLUL","MLeNN","MLTL"),])) { times[i,2] <- as.numeric(times[i,2]) + timeNeighbors }
        for (i in rownames(times[times$algorithm %in% c("MLSMOTE","MLSOL","MLUL","MLeNN","MLTL"),])) { times[i,2] <- as.numeric(times[i,2]) + timeTableVDM }

        message(paste("End of execution. Generated MLDs stored under directory",outputDirectory))

        times

      }

    } else {

      tableVDM <- NULL
      timeTableVDM <- 0
      neighbors <- NULL
      neighbors2 <- NULL
      timeNeighbors <- 0

      message(paste("Calculating structures for dataset", D$name, ", if necessary. Once this is done, algorithms will be applied faster"))

      if (length(intersect(algorithms, c("MLSMOTE","MLeNN","MLTL","MLSOL","MLUL"))) > 0) {

        message(paste("Calculating VDM table for dataset", D$name))
        startTime <- Sys.time()
        tableVDM <- calculateTableVDM(D)
        endTime <- Sys.time()
        timeTableVDM <- as.numeric(endTime - startTime, units="secs")
        message(paste("Time taken (in seconds):",timeTableVDM))

      }

      if (length(intersect(algorithms, c("MLeNN", "MLTL","MLSOL","MLUL"))) > 0) {

        message(paste("Calculating neighbors structure for dataset", D$name, ". Once this is done, algorithms will be applied faster"))
        startTime <- Sys.time()
        neighbors <- getAllNeighbors(D, c(1:D$measures$num.instances)[D$dataset$.labelcount > 0], tableVDM)
        neighbors2 <- getAllNeighbors2(neighbors, unique(unlist(mldr.resampling.env$.mldrApplyFun1(D$labels[D$labels$IRLbl < D$measures$meanIR,]$index, function(x) { c(1:D$measures$num.instances)[D$dataset[x]==1] }, mc.cores=mldr.resampling.env$.numCores))), k)
        neighbors <- mldr.resampling.env$.mldrApplyFun1(neighbors, function(x) { x[1:k+1] }, mc.cores=mldr.resampling.env$.numCores)
        endTime <- Sys.time()
        timeNeighbors <- as.numeric(endTime - startTime, units="secs")
        message(paste("Time taken (in seconds):",timeNeighbors))

      }

      for (a in algorithms) {

        time <- executeAlgorithm(D, a, P, k, TH, outputDirectory, neighbors, neighbors2, tableVDM)

        if (!is.null(time)) {
          times[nrow(times) + 1,] <- c(a,time)
        }

      }

      #Add structure generation times
      for (i in rownames(times[times$algorithm %in% c("MLSOL","MLUL","MLeNN","MLTL"),])) { times[i,2] <- as.numeric(times[i,2]) + timeNeighbors }
      for (i in rownames(times[times$algorithm %in% c("MLSMOTE","MLSOL","MLUL","MLeNN","MLTL"),])) { times[i,2] <- as.numeric(times[i,2]) + timeTableVDM }

      message(paste("End of execution. Generated MLDs stored under directory",outputDirectory))

      times

    }

  }

}



#' Set the number of cores available for parallel computing
#'
#' @param n The new value for the number of cores
#'
#' @return No return value, called in order to change the number of cores
#'
#' @examples
#' \donttest{
#' setNumCores(8)
#' }
#' @export
setNumCores <- function(n) {
  assign('mldr.resampling.env$.numCores', n, mldr.resampling.env)
}



#' Get the number of cores available for parallel computing
#'
#' @return The number of cores available for parallel computing
#'
#' @examples
#' getNumCores()
#'
#' @export
getNumCores <- function() {
  mldr.resampling.env$.numCores
}



#' Enable/Disable parallel computing
#'
#' @param beParallel A boolean indicating if parallel computing is to be enabled (TRUE) or disabled (FALSE)
#'
#' @return No return value, called in order to enable parallel computing
#'
#' @examples
#' setParallel(TRUE)
#'
#' @export
setParallel <- function(beParallel) {
  if (!beParallel) {
    assign('mldr.resampling.env$.mldrApplyFun2', function(x, l, mc.cores) { pbapply::pblapply(x,l) }, mldr.resampling.env)
    message("Parallel computing disabled")
  } else {
    if (requireNamespace("parallel", quietly = TRUE)) {
      setNumCores(parallel::detectCores())
      assign('mldr.resampling.env$.mldrApplyFun2', parallel::mclapply, mldr.resampling.env)
      message(paste("Parallel computing enabled on all",getNumCores(),"available cores. Use function setNumCores if you wish to modify it"))
    } else {
      message("You have to install package parallel in order to enable parallel computing")
    }
  }
}
