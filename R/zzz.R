.onAttach <- function(...) {
  packageStartupMessage("Enter setParallel(TRUE) to enable parallel computing")
}


.onLoad <- function(...) {
  .numCores <<- 1
  .mldrApplyFun1 <<- function(x, l, mc.cores) { lapply(x,l) }
  .mldrApplyFun2 <<- function(x, l, mc.cores) { pbapply::pblapply(x,l) }
}


.onUnload <- function(...) {
  rm(.mldrApplyFun1)
  rm(.mldrApplyFun2)
  rm(.numCores)
}
