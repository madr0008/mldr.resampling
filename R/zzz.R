mldr.resampling.env <- new.env(parent=emptyenv())


.onAttach <- function(...) {
  packageStartupMessage("Enter setParallel(TRUE) to enable parallel computing")
}


.onLoad <- function(...) {

  assign('.numCores', 1, mldr.resampling.env)
  assign('.mldrApplyFun1', function(x, l, mc.cores) { lapply(x,l) }, mldr.resampling.env)
  assign('.mldrApplyFun2', function(parV, parF, mc.cores, parL, parEnv) { pbapply::pblapply(parV, parF) }, mldr.resampling.env)

}


.onUnload <- function(...) {
  rm(mldr.resampling.env)
}
