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
#'
#' @return A synthetic sample derived from the one passed as a parameter and its neighbors
#' @examples
#' \dontrun{
#' library(mldr)
#' initTypes(C,c(34,65,121),3,bibtex)
#' }
initTypes <- function(C, neighbors, k, minoritary, D) {

  t <- lapply(as.numeric(rownames(D$dataset)), function(i) {
    unlist(lapply(D$labels$index, function(j) {
      if (D$dataset[i,j] != minoritary[j - D$measures$num.inputs]) {
        5
      } else {
        if (C[[as.character(i)]][j - D$measures$num.inputs] < 1/3) {
          1
        } else if (C[[as.character(i)]][j - D$measures$num.inputs] < 2/3) {
          2
        } else if (C[[as.character(i)]][j - D$measures$num.inputs] < 1) {
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
    for (i in 1:D$measures$num.instances) {
      for (j in 1:D$measures$num.labels) {
        if (t[[i]][j] == 3) {
          for (m in neighbors[[as.character(i)]]) {
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

  attributes <- as.numeric(unlist(lapply(D$attributesIndexes, function(i) { #Attributes
    ifelse(D$attributes[[i]] %in% c("numeric", "Date"),
           D$dataset[s,i] + (D$dataset[r,i] - D$dataset[s,i])*stats::runif(1, 0, 1), #Numeric attributes
           sample(c(D$dataset[s,i], D$dataset[r,i]), size = 1)) #Non numeric attributes
  })))

  #Calculate distances between attributes
  d_s <- sum((attributes - as.numeric(unlist(D$dataset[s,D$attributesIndexes])))^2)
  d_r <- sum((attributes - as.numeric(unlist(D$dataset[r,D$attributesIndexes])))^2)

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

  stats::setNames(lapply(d, function(i) {
    activeLabels <- D$labels[which(D$dataset[i,D$labels$index] %in% 1),1]
    if (length(activeLabels)>0) {
      getNN(i, d, ifelse(length(activeLabels)==1,activeLabels,sample(activeLabels,1)), D, k)
    }
  }), d)

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

  stats::setNames(lapply(d, function(i) {
    unlist(lapply(D$labels$index, function(j) {
      (1/k) * sum(unlist(lapply(neighbors[[as.character(i)]], function(m) {
        ifelse(D$dataset[i,j]==D$dataset[m,j],0,1)
      })))
    }))
  }), d)

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

  stats::setNames(lapply(d, function(i) {
    unlist(lapply(D$labels$index, function(j) {
      if ((D$dataset[i,j]==minoritary[j - D$measures$num.inputs]) & (C[[as.character(i)]][j - D$measures$num.inputs]<1)) {
        numerator <- C[[as.character(i)]][j - D$measures$num.inputs]
        denominator <- sum(unlist(lapply(d, function(x) {
          ifelse(D$dataset[x,j]==minoritary[j - D$measures$num.inputs],C[[as.character(x)]][j - D$measures$num.inputs],0)*ifelse(C[[as.character(x)]][j - D$measures$num.inputs]<1,1,0)
        })))
        numerator/denominator
      } else {
        -1
      }
    }))
  }), d)

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

  stats::setNames(lapply(d, function(i) {
    sum(S[[as.character(i)]][!S[[as.character(i)]] %in% -1])
  }), d)

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

  stats::setNames(lapply(d, function(i) {
    d[ceiling(which(unlist(neighbors)==i)/k)]
  }), d)

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

  stats::setNames(lapply(d, function(i) {

    sum(unlist(lapply(D$labels$index, function(j) {

        ifelse(length(rNeighbors[[as.character(i)]]) > 0,

        sum(unlist(lapply(rNeighbors[[as.character(i)]], function(m) {

          ifelse(S[[as.character(m)]][j - D$measures$num.inputs]==-1,0,ifelse(D$dataset[i,j]==D$dataset[m,j],1,-1)*S[[as.character(m)]][j - D$measures$num.inputs])

        }))) / length(rNeighbors[[as.character(i)]]), 0)

    })))

  }), d)

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

  v <- lapply(d, function(i) {
         w[[as.character(i)]] + u[[as.character(i)]]
       })

  minimum <- min(unlist(v))

  stats::setNames(lapply(v, function(x) {
    x - minimum
  }), d)

}



#' Auxiliary function used by resample. It executes an algorithm, given as a string, and stores the resulting MLD in a arff file
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param a String with the name of the algorithm to be applied.
#' @param P Percentage in which the original dataset is increased/decreased (if required by the algorithm)
#' @param k Number of neighbors taken into account for each instance (if required by the algorithm)
#' @param TH Threshold for the Hamming Distance in order to consider an instance different to another one (if required by the algorithm)
#' @param outputDirectory Route with the directory where the generated ARFF file will be stored
#'
#' @examples
#' \dontrun{
#' library(mldr)
#' executeAlgorithm(bibtex, "MLSMOTE", k=3)
#' }
executeAlgorithm <- function(D, a, P, k, TH, outputDirectory) {

  if (!(a %in% c("LPROS", "LPRUS", "MLROS", "MLRUS", "MLRkNNOS", "MLSMOTE", "MLSOL", "MLUL", "MLeNN", "MLTL", "REMEDIAL"))) {
    print(paste("Error: There is no algorithm named", a))
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
      d <- f(D, P, k)
      endTime <- Sys.time()
    } else if (a == "MLeNN") {
      name <- paste(D$name, a, "TH", TH, "k", k, sep = "_")
      print(paste("Running",a,"on",D$name,"with TH =",TH,"and k =",k))
      startTime <- Sys.time()
      d <- f(D, TH, k)
      endTime <- Sys.time()
    } else if (a == "MLTL") {
      name <- paste(D$name, a, "TH", TH, sep = "_")
      print(paste("Running",a,"on",D$name,"with TH =",TH))
      startTime <- Sys.time()
      d <- f(D, TH)
      endTime <- Sys.time()
    } else { #REMEDIAL
      name <- a
      print(paste("Running",a,"on",D$name))
      startTime <- Sys.time()
      d <- f(D)
      endTime <- Sys.time()
    }

    print(paste("Time taken (in seconds):",as.numeric(endTime - startTime, units="secs")))

    mldr::write_arff(d, paste(outputDirectory, name, sep="/"))

  }

}



#' Interface function of the package. It executes one or several, given as a strings, and stores the resulting MLDs in arff files
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param algorithms String, or string vector, with the name(s) of the algorithm(s) to be applied.
#' @param P Percentage in which the original dataset is increased/decreased, if required by the algorithm(s). Defaults to 25
#' @param k Number of neighbors taken into account for each instance, if required by the algorithm(s). Defaults to 3
#' @param TH Threshold for the Hamming Distance in order to consider an instance different to another one, if required by the algorithm(s). Defaults to 0.5
#' @param params Dataframe with 4 columns: name of the algorithm, P, k and TH, in order to execute several algorithms with different values for their parameters
#' @param outputDirectory Route with the directory where generated ARFF files will be stored. Defaults to the working directory
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

  if (missing(D)) {

    print("Please, provide a mld object as the original dataset")

  } else {

    if (missing(algorithms)) {

      if (missing(params)) {

        print("Please, specify the algorithms to be applied, either with the parameter algorithms or with params")

      } else {

        for(i in 1:nrow(params)) {

          executeAlgorithm(D,params[i,1],params[i,2],params[i,3],params[i,4])

        }

        print(paste("End of execution. Generated MLDs stored under directory",outputDirectory))

      }

    } else {

      for (a in algorithms) {

        executeAlgorithm(D, a, P, k, TH, outputDirectory)

      }

      print(paste("End of execution. Generated MLDs stored under directory",outputDirectory))

    }

  }

}
