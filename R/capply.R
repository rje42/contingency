##' Apply function over tables
##' 
##' Apply a function to each contingency table in a `tables` object.
##' 
##' @param x object of class `tables`
##' @param f function to apply to each table
##' @param ... additional arguments to `f`
##' 
##' @return a vector, matrix or list of outputs from the function `f`.
##' 
##' @export
capply <- function(x, f, ...) {
  if (is_rev(x)) {
    x2 <- as.array(x)
    tmp <- apply(x2, length(dims(x2)), f, ...)
  }
  else tmp <- apply(as.array(x), 1, f, ...)
  if (is.matrix(tmp)) return(t(tmp))
  else tmp
}

##' Turn distributions into tables
##' 
##' @param n number of distributions to generate
##' @param f function that generates a probability distribution
##' @param ... arguments to `f`
##' @param rev logical: should output move through each table fastest?
##' 
##' @return a tables object containing the outputs of `f`
##' 
##' @export
repTables <- function (n, f, ..., rev=FALSE) {
  args <- list(...)
  out <- replicate(n, do.call(f, args))
  if (is.null(dim(out))) stop("f should give distributions")
  
  if (rev) {
    do <- dim(out)
    out <- as_tables(out, tdim=do[-length(do)])
  }
  else {
    out <- aperm(out, c(length(dim(out)), seq_along(dim(out)[-1])))
    out <- as_tables(out, tdim=dim(out)[-1])
  }
  
  return(out)
}
