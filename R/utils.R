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

  unlist(mldrApplyFun1(rest, function(y) {
    ifelse(y == sample,
           Inf, #In order not to choose its own
           sqrt( #Square root because euclidean distance
             sum( #Summing distances between numeric and non numeric attributes
               sum( #For numeric attributes: square of the difference
                 (D$dataset[sample,D$attributesIndexes[D$attributes[1:D$measures$num.inputs]=="numeric"]] - D$dataset[y,D$attributesIndexes[D$attributes[1:D$measures$num.inputs]=="numeric"]])^2
               ),
               ifelse(sum(D$attributes[1:D$measures$num.inputs]!="numeric") > 0,
                 sum( #For non numeric attributes: Value Difference Measure (VDM)
                   unlist(mldrApplyFun1(D$attributesIndexes[D$attributes[1:D$measures$num.inputs]!="numeric"], function(x) {
                     table1 <- table((D$dataset[D$dataset[x] == D$dataset[sample,x],])[label])/(table(D$dataset[x])[[D$dataset[sample,x]]])
                     table2 <- table((D$dataset[D$dataset[x] == D$dataset[y,x],])[label])/(table(D$dataset[x])[[D$dataset[y,x]]])
                     sum(
                       abs(ifelse(length(table1 == 1), ifelse(names(table1) == "0", stats::setNames(c(table1, 0), c("0","1")), stats::setNames(c(0, table1), c("0","1"))), table1) - ifelse(length(table2 == 1), ifelse(names(table2) == "0", stats::setNames(c(table2, 0), c("0","1")), stats::setNames(c(0, table2), c("0","1"))), table2))
                     )
                   }))
                 ), 0)
             )
           )
    )
  }))

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
    mldrApplyFun1(D$attributesIndexes[1:D$measures$num.inputs], function(i) { #Attributes
      ifelse(D$attributes[[i]] %in% c("numeric", "Date"),
             D$dataset[seedInstance,i] + (D$dataset[refNeigh,i] - D$dataset[seedInstance,i])*stats::runif(1, 0, 1), #Numeric attributes
             utils::tail(names(sort(table(D$dataset[neighbors, i]))), 1)) #Non numeric attributes
    }),
    unlist(mldrApplyFun1(mldrApplyFun1(D$dataset[c(seedInstance, neighbors),D$labels$index], sum), function(x) { #Labels
      ifelse(x > ((length(neighbors)+1)/2), 1, 0)
    }),
    rep(NA, length(D$dataset) - D$measures$num.attributes) #Other measures like labelcount, SCUMBLE
  ))

}



#' Auxiliary function used by MLeNN. Computes the Hamming Distance between two instances
#'
#' @param x Index of sample 1
#' @param y Index of sample 2
#' @param D mld \code{mldr} object in which the instances are located
#'
#' @return The Hamming Distance between the instances
#' @examples
#' \dontrun{
#' library(mldr)
#' adjustedHammingDist(1,2,bibtex)
#' }
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
#' @examples
#' \dontrun{
#' library(mldr)
#' initTypes(C,c(34,65,121),3,bibtex)
#' }
initTypes <- function(C, neighbors, k, minoritary, D, d) {

  t <- mldrApplyFun2(d, function(i) {
    unlist(mldrApplyFun1(D$labels$index, function(j) {
      if (D$dataset[i,j] != minoritary[j - D$measures$num.inputs]) {
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
    }))
  })

  change <- TRUE
  while (change) {
    change <- FALSE
    for (i in d) {
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
#' @examples
#' \dontrun{
#' library(mldr)
#' generateInstanceMLSOL(25,28,c(1,3),bibtex)
#' }
generateInstanceMLSOL <- function (seedInstance, refNeigh, t, D) {

  s <- seedInstance
  r <- refNeigh

  attributes <- as.numeric(unlist(mldrApplyFun1(D$attributesIndexes[1:D$measures$num.inputs], function(i) { #Attributes
    ifelse(D$attributes[[i]] %in% c("numeric", "Date"),
           D$dataset[s,i] + (D$dataset[r,i] - D$dataset[s,i])*stats::runif(1, 0, 1), #Numeric attributes
           sample(c(D$dataset[s,i], D$dataset[r,i]), size = 1)) #Non numeric attributes
  })))

  #Calculate distances between attributes
  d_s <- sum((attributes - as.numeric(unlist(D$dataset[s,D$attributesIndexes[1:D$measures$num.inputs]])))^2)
  d_r <- sum((attributes - as.numeric(unlist(D$dataset[r,D$attributesIndexes[1:D$measures$num.inputs]])))^2)

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

  c(attributes, labels, rep(NA, length(D$dataset) - D$measures$num.attributes))

}



#' Auxiliary function used by MLSOL and MLUL. Computes the kNN of every instance in a dataset
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param k Number of neighbors to be included for every instance
#'
#' @return A list of vectors with the indexes of the neighbors for each instance
#' @examples
#' \dontrun{
#' library(mldr)
#' getAllNeighbors(bibtex, 3)
#' }
getAllNeighbors <- function(D, d, k) {

  mldrApplyFun2(d, function(i) {
    activeLabels <- D$labels[which(D$dataset[i,D$labels$index] %in% 1),1]
    if (length(activeLabels) > 0) {
      getNN(i, d, ifelse(length(activeLabels)==1,activeLabels,sample(activeLabels,1)), D, k)
    }
  })

}



#' Auxiliary function used by MLSOL and MLUL. For each instance in the dataset, we compute, for each label, the proportion of neighbors having an opposite class with respect to the proper instance
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param neighbors Structure with the neighbors of every instance in the dataset
#' @param k Number of neighbors taken into account for each instance
#'
#' @return A structure with the proportion of neighbors having an opposite class with respect to an instance and label
#' @examples
#' \dontrun{
#' library(mldr)
#' getC(bibtex, rownames(bibtex), neighbors, 3)
#' }
getC <- function(D, d, neighbors, k) {

  mldrApplyFun2(d, function(i) {
    unlist(mldrApplyFun1(D$labels$index, function(j) {
      (1/k) * sum(unlist(mldrApplyFun1(neighbors[i], function(m) {
        ifelse(D$dataset[i,j]==D$dataset[m,j],0,1)
      })))
    }))
  })

}



#' Auxiliary function used by MLSOL and MLUL. For non outlier instances, it aggregates the values of C, taking into account the global class imbalance
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param C Structure with the proportion of neighbors having an opposite class with respect to an instance and label
#' @param minoritary Vector with the minoritary class of each label (normally, 1)
#'
#' @return A structure with the proportion of neighbors having an opposite class with respect to an instance and label, normalized by the global class imbalance
#' @examples
#' \dontrun{
#' library(mldr)
#' getS(bibtex, rownames(bibtex), C, rep(1, bibtex$measures$num.labels))
#' }
getS <- function(D, d, C, minoritary) {

  mldrApplyFun2(d, function(i) {
    unlist(mldrApplyFun1(D$labels$index, function(j) {
      if ((D$dataset[i,j]==minoritary[j - D$measures$num.inputs]) & (C[[i]][j - D$measures$num.inputs]<1)) {
        numerator <- C[[i]][j - D$measures$num.inputs]
        denominator <- sum(unlist(mldrApplyFun1(d, function(x) {
          ifelse(D$dataset[x,j]==minoritary[j - D$measures$num.inputs],C[[x]][j - D$measures$num.inputs],0)*ifelse(C[[x]][j - D$measures$num.inputs]<1,1,0)
        })))
        numerator/denominator
      } else {
        -1
      }
    }))
  })

}



#' Auxiliary function used by MLSOL and MLUL. For non outlier instances, it aggregates the values of S for each label
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param S Structure with the proportion of neighbors having an opposite class with respect to an instance and label, normalized by the global class imbalance
#'
#' @return A vector of weights to be considered when oversampling for each instance
#' @examples
#' \dontrun{
#' library(mldr)
#' getW(bibtex, rownames(bibtex), S)
#' }
getW <- function(D, d, S) {

  mldrApplyFun2(d, function(i) {
    sum(S[[i]][!S[[i]] %in% -1])
  })

}



#' Auxiliary function used by MLUL. For each instance in the dataset, given the neighbors structure, we compute its reverse nearest neighbors
#'
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param neighbors Structure with the neighbors of every instance in the dataset
#' @param k Number of neighbors to be considered
#'
#' @return A list of vectors with the indexes of the reverse nearest neighbors of every instance in the dataset
#' @examples
#' \dontrun{
#' library(mldr)
#' getAllReverseNeighbors(rownames(bibtex), neighbors, 3)
#' }
getAllReverseNeighbors <- function(d, neighbors, k) {

  mldrApplyFun2(d, function(i) {
    d[ceiling(which(unlist(neighbors)==i)/k)]
  })

}



#' Auxiliary function used by MLUL. It computes the influence of each instance with respect to its reverse neighbors
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param rNeighbors Structure with the reverse nearest neighbors of each instance of the dataset
#' @param S Structure with the proportion of neighbors having an opposite class with respect to an instance and label, normalized by the global class imbalance
#'
#' @return A list of values of influence for each instance with respect to its reverse neighbors
#' @examples
#' \dontrun{
#' library(mldr)
#' getU(bibtex, rownames(bibtex), rNeighbors, S)
#' }
getU <- function(D, d, rNeighbors, S) {

  mldrApplyFun2(d, function(i) {

    sum(unlist(mldrApplyFun1(D$labels$index, function(j) {

        ifelse(length(rNeighbors[[i]]) > 0,

        sum(unlist(mldrApplyFun1(rNeighbors[[i]], function(m) {

          ifelse(S[[m]][j - D$measures$num.inputs]==-1,0,ifelse(D$dataset[i,j]==D$dataset[m,j],1,-1)*S[[m]][j - D$measures$num.inputs])

        }))) / length(rNeighbors[i]), 0)

    })))

  })

}



#' Auxiliary function used by MLUL. It calculates, for each instance, how important it is in the dataset
#'
#' @param d Vector with the instances of the dataset which have one or more label active (ideally, all of them)
#' @param w List of weights for each instance
#' @param u List of influences in reverse neighbors for each instance
#'
#' @return A list with the values of importance of each instance in the dataset
#' @examples
#' \dontrun{
#' library(mldr)
#' getV(rownames(bibtex), w, u)
#' }
getV <- function(d, w, u) {

  v <- mldrApplyFun2(d, function(i) {
         w[[i]] + u[[i]]
       })

  minimum <- min(unlist(v))

  mldrApplyFun1(v, function(x) {
    x - minimum
  })

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
#'
#' @return Time (in seconds) taken to execute the algorithm (NULL if no algorithm was executed)
#'
#' @examples
#' \dontrun{
#' library(mldr)
#' executeAlgorithm(bibtex, "MLSMOTE", k=3)
#' }
executeAlgorithm <- function(D, a, P, k, TH, outputDirectory, neighbors, neighbors2) {

  if (!(a %in% c("LPROS", "LPRUS", "MLROS", "MLRUS", "MLRkNNOS", "MLSMOTE", "MLSOL", "MLUL", "MLeNN", "MLTL", "REMEDIAL"))) {
    print(paste("Error: There is no algorithm named", a))
    NULL
  } else {

    f <- get(a)
    if (a %in% c("LPROS", "LPRUS", "MLROS", "MLRUS")) {
      name <- paste(D$name, a, "P", P, sep = "_")
      print(paste("Running",a,"on",D$name,"with P =",P))
      startTime <- Sys.time()
      d <- f(D, P)
      endTime <- Sys.time()
    } else if (a %in% c("MLRkNNOS", "MLSMOTE")) {
      name <- paste(D$name, a, "k", k, sep = "_")
      print(paste("Running",a,"on",D$name,"with k =",k))
      startTime <- Sys.time()
      d <- f(D, k)
      endTime <- Sys.time()
    } else if (a %in% c("MLSOL", "MLUL")) {
      name <- paste(D$name, a, "P", P, "k", k, sep = "_")
      print(paste("Running",a,"on",D$name,"with P =",P,"and k =",k))
      startTime <- Sys.time()
      d <- f(D, P, k, neighbors)
      endTime <- Sys.time()
    } else if (a == "MLeNN") {
      name <- paste(D$name, a, "TH", TH, "k", k, sep = "_")
      print(paste("Running",a,"on",D$name,"with TH =",TH,"and k =",k))
      startTime <- Sys.time()
      d <- f(D, TH, k, neighbors2)
      endTime <- Sys.time()
    } else if (a == "MLTL") {
      name <- paste(D$name, a, "TH", TH, sep = "_")
      print(paste("Running",a,"on",D$name,"with TH =",TH))
      startTime <- Sys.time()
      d <- f(D, TH, neighbors2)
      endTime <- Sys.time()
    } else { #REMEDIAL
      name <- a
      print(paste("Running",a,"on",D$name))
      startTime <- Sys.time()
      d <- f(D)
      endTime <- Sys.time()
    }

    time <- as.numeric(endTime - startTime, units="secs")

    print(paste("Time taken (in seconds):",time))

    mldr::write_arff(d, paste(outputDirectory, name, sep="/"))

    time

  }

}



#' Interface function of the package. It executes one or several, given as a strings, and stores the resulting MLDs in arff files
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param algorithms String, or string vector, with the name(s) of the algorithm(s) to be applied.
#' @param P Percentage in which the original dataset is increased/decreased, if required by the algorithm(s). Defaults to 25
#' @param k Number of neighbors taken into account for each instance, if required by the algorithm(s). Defaults to 3
#' @param TH Threshold for the Hamming Distance in order to consider an instance different to another one, if required by the algorithm(s). Defaults to 0.5
#' @param params Dataframe with 4 columns: name of the algorithm, P, k and TH, in that order, to execute several algorithms with different values for their parameters
#' @param outputDirectory Route with the directory where generated ARFF files will be stored. Defaults to the working directory
#'
#' @return Dataframe with times (in seconds) taken in to execute each algorithm
#'
#' @examples
#' \dontrun{
#' library(mldr)
#' resample(bibtex, "MLSMOTE", k=3)
#' resample(bibtex, c("MLSOL", "MLeNN"), P=30, k=5, TH=0.4)
#' resample(bibtex, params)
#' }
#' @export
resample <- function(D, algorithms, P=25, k=3, TH=0.5, params, outputDirectory=getwd()) {

  times <- data.frame(matrix(nrow = 0, ncol = 2))
  colnames(times) <- c("algorithm", "time")

  if (missing(D)) {

    print("Please, provide a mld object as the original dataset")
    NULL

  } else {

    if (missing(algorithms)) {

      if (missing(params)) {

        print("Please, specify the dataset and algorithms to be applied, either with the parameter algorithms or with params")
        NULL

      } else {

        neighbors <- NULL
        timeNeighbors <- 0
        neighbors2 <- NULL
        timeNeighbors2 <- 0

        if (("MLSOL" %in% unlist(params[1]) | "MLUL" %in% unlist(params[1])) & (sum(ifelse(is.na(table(params[1])["MLSOL"]),0,table(params[1])["MLSOL"]), ifelse(is.na(table(params[1])["MLUL"]),0,table(params[1])["MLUL"])) > 1)) {

          print(paste("Calculating neighbors structures for dataset", D$name, ". Once this is done, all the algorithms will be applied faster"))
          startTime <- Sys.time()
          neighbors <- getAllNeighbors(D, c(1:D$measures$num.instances)[D$dataset$.labelcount > 0], max(params[params[,1] %in% c("MLSOL","MLUL"),3]))
          endTime <- Sys.time()
          timeNeighbors <- as.numeric(endTime - startTime, units="secs")
          print(paste("Time taken (in seconds):",timeNeighbors))

        }

        if (("MLeNN" %in% unlist(params[1]) | "MLTL" %in% unlist(params[1])) & (sum(ifelse(is.na(table(params[1])["MLeNN"]),0,table(params[1])["MLeNN"]), ifelse(is.na(table(params[1])["MLTL"]),0,table(params[1])["MLTL"])) > 1)) {

          print(paste("Calculating second neighbors structure for dataset", D$name, ". Once this is done, algorithms MLeNN and MLTL will be applied faster"))
          startTime <- Sys.time()
          neighbors2 <- getAllNeighbors(D, unique(unlist(mldrApplyFun1(D$labels[D$labels$IRLbl < D$measures$meanIR,]$index, function(x) c(1:D$measures$num.instances)[D$dataset[x]==1]))), max(params[params[,1] %in% c("MLSOL","MLUL"),3]))
          endTime <- Sys.time()
          timeNeighbors2 <- as.numeric(endTime - startTime, units="secs")
          print(paste("Time taken (in seconds):",timeNeighbors2))

        }

        for(i in 1:nrow(params)) {

          time <- executeAlgorithm(D,params[i,1],params[i,2],params[i,3],params[i,4], neighbors, neighbors2)

          if (!is.null(time)) {
            times[nrow(times) + 1,] <- c(a,time)
          }

        }

        #Add neighbors structure generation time
        for (i in rownames(times[times$algorithm %in% c("MLSOL","MLUL"),])) { times[i,2] <- as.numeric(times[i,2]) + timeNeighbors }
        for (i in rownames(times[times$algorithm %in% c("MLeNN","MLTL"),])) { times[i,2] <- as.numeric(times[i,2]) + timeNeighbors2 }

        print(paste("End of execution. Generated MLDs stored under directory",outputDirectory))

        times

      }

    } else {

      neighbors <- NULL
      timeNeighbors <- 0
      neighbors2 <- NULL
      timeNeighbors2 <- 0

      if (("MLSOL" %in% algorithms | "MLUL" %in% algorithms) & (sum(ifelse(is.na(table(algorithms)["MLSOL"]),0,table(algorithms)["MLSOL"]), ifelse(is.na(table(algorithms)["MLUL"]),0,table(algorithms)["MLUL"])) > 1)) {

        print(paste("Calculating neighbors structure for dataset", D$name, ". Once this is done, algorithms MLSOL and MLUL will be applied faster"))
        startTime <- Sys.time()
        neighbors <- getAllNeighbors(D, c(1:D$measures$num.instances)[D$dataset$.labelcount > 0], k)
        endTime <- Sys.time()
        timeNeighbors <- as.numeric(endTime - startTime, units="secs")
        print(paste("Time taken (in seconds):",timeNeighbors))

      }

      if (("MLeNN" %in% algorithms | "MLTL" %in% algorithms) & (sum(ifelse(is.na(table(algorithms)["MLeNN"]),0,table(algorithms)["MLeNN"]), ifelse(is.na(table(algorithms)["MLTL"]),0,table(algorithms)["MLTL"])) > 1)) {

        print(paste("Calculating second neighbors structure for dataset", D$name, ". Once this is done, algorithms MLeNN and MLTL will be applied faster"))
        startTime <- Sys.time()
        neighbors2 <- getAllNeighbors(D, unique(unlist(mldrApplyFun1(D$labels[D$labels$IRLbl < D$measures$meanIR,]$index, function(x) c(1:D$measures$num.instances)[D$dataset[x]==1]))), k)
        endTime <- Sys.time()
        timeNeighbors2 <- as.numeric(endTime - startTime, units="secs")
        print(paste("Time taken (in seconds):",timeNeighbors2))

      }

      for (a in algorithms) {

        time <- executeAlgorithm(D, a, P, k, TH, outputDirectory, neighbors, neighbors2)

        if (!is.null(time)) {
          times[nrow(times) + 1,] <- c(a,time)
        }

      }

      #Add neighbors structure generation time
      for (i in rownames(times[times$algorithm %in% c("MLSOL","MLUL"),])) { times[i,2] <- as.numeric(times[i,2]) + timeNeighbors }
      for (i in rownames(times[times$algorithm %in% c("MLeNN","MLTL"),])) { times[i,2] <- as.numeric(times[i,2]) + timeNeighbors2 }

      print(paste("End of execution. Generated MLDs stored under directory",outputDirectory))

      times

    }

  }

}
