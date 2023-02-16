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



#' @title Multi-label oversampling based on local label imbalance (MLSOL)
#'
#' @description This function implements the MLSOL algorithm. It is a preprocessing algorithm for imbalanced multilabel datasets,
#' whose aim is to identify instances with minoritary labels, and generate synthetic instances based on their neighbor instances.
#'
#' @source Liu, B., Blekas, K., & Tsoumakas, G. (2022). Multi-label sampling based on local label imbalance. Pattern Recognition, 122, 108294.
#'
#' @param D mld \code{mldr} object with the multilabel dataset to preprocess
#' @param P Percentage in which the original dataset is increased
#' @param k Number of neighbors to be considered when computing the neighbors of an instance
#'
#' @return A mld object containing the preprocessed multilabel dataset
#' @examples
#' \dontrun{
#' library(mldr)
#' MLSOL(bibtex, 3)
#' }
#' @export
MLSOL <- function(D, P, k) {

  minoritary <- unlist(lapply(D$labels$freq, function(x) ifelse(x<0.5,1,0)))

  neighbors <- stats::setNames(lapply(as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,])), function(i) {
    activeLabels <- D$labels[which(D$dataset[i,D$labels$index] %in% 1),1]
    if (length(activeLabels)>0) {
      order(calculateDistances(i, as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,])), ifelse(length(activeLabels)==1,activeLabels,sample(activeLabels,1)), D))[1:k+1]
    }
  }), as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,])))


  C <- stats::setNames(lapply(as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,])), function(i) {
    unlist(lapply(D$labels$index, function(j) {
      (1/k) * sum(unlist(lapply(neighbors[[as.character(i)]], function(m) {
        ifelse(D$dataset[i,j]==D$dataset[m,j],0,1)
      })))
    }))
  }), as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,])))

  S <- stats::setNames(lapply(as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,])), function(i) {
    unlist(lapply(D$labels$index, function(j) {
      if ((D$dataset[i,j]==minoritary[j - D$measures$num.inputs]) & (C[[as.character(i)]][j - D$measures$num.inputs]<1)) {
        numerator <- C[[as.character(i)]][j - D$measures$num.inputs]
        denominator <- sum(unlist(lapply(as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,])), function(x) {
          ifelse(D$dataset[x,j]==minoritary[j - D$measures$num.inputs],C[[as.character(x)]][j - D$measures$num.inputs],0)*ifelse(C[[as.character(x)]][j - D$measures$num.inputs]<1,1,0)
        })))
        numerator/denominator
      } else {
        -1
      }
    }))
  }), as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,])))

  w <- stats::setNames(lapply(as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,])), function(i) {
    sum(S[[as.character(i)]][!S[[as.character(i)]] %in% -1])
  }), as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,])))

  t <- initTypes(C, neighbors, k, minoritary)

  genNum <- (D$measures$num.instances/100) * P
  seedInstances <- sample(as.numeric(rownames(D$dataset[D$dataset$.labelcount > 0,])), size=genNum, replace=TRUE, prob=w)
  newSamples <- lapply(seedInstances, function(i) {
    refNeigh <- sample(neighbors[[as.character(i)]], size=1)
    generateInstanceMLSOL(i, refNeigh, t, D)
  })

  mldr::mldr_from_dataframe(rbind(D$dataset, lapply(stats::setNames(as.data.frame(do.call(rbind, newSamples[-1])), names(D$dataset)), unlist)), D$labels$index, D$attributes, D$name)

}
