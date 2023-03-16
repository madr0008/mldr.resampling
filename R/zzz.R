.onLoad <- function(...) {
  if (requireNamespace("parallel", quietly = TRUE)) {
    mldrApplyFun1 <<- parallel::mclapply
    mldrApplyFun2 <<- parallel::mclapply
  } else {
    mldrApplyFun1 <<- lapply
    mldrApplyFun2 <<- pbapply::pblapply
  }
}

.onUnload <- function(...) {
  rm(mldrApplyFun1)
  rm(mldrApplyFun2)
}
