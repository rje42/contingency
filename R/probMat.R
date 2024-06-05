##' Generate matrix of (conditional) probability distributions
##' 
##' @description Generates discrete probability distributions in a matrix.
##' 
##' @param n number of distributions
##' @param dim dimension of contingency table for distributions
##' @param d number of dimensions of table
##' @param alpha parameter to use in dirichlet distribution
##' @param condition which dimensions should be conditioned upon
##' 
##' @details Returns an object of class `tables` consisting of 
##' discrete probability distributions.  Each distribution is assumed to be a 
##' contingency table of dimension `dim`, and the probabilities
##' are generated using a Dirichlet distribution with parameters all 
##' equal to `alpha`.
##' 
##' @return A `tables` object containing random distributions.
##' 
##' @examples
##' dat <- rprobMat(10, c(2,2,2))
##' 
##' @export 
rprobMat <- function(n, dim, d, alpha=1) {
  ## if dimension vector shorter than length d, recycle (with warning
  ## if necessary)
  if (missing(d)) d <- length(dim)
  else if (length(dim) < d) dim <- dim * rep.int(1L, d)
  else if (length(dim) > d) stop("More than 'd' dimensions supplied")
  
  if (any(dim < 0)) stop("Dimensions must be non-negative")
  if (any(alpha < 0)) stop("Parameters must be non-negative")
  
  k <- prod(dim)
  out <- matrix(rgamma(n*k, alpha, 1), n, k)
#  out <- out/rep(.colSums(out, k, n), each=k)
  out <- out/.rowSums(out, n, k)
  
  class(out) <- "tables"
  attr(out, "tdim") <- dim
  attr(out, "conditional") <- integer(0)
#  dim(out) <- c(dim, n)
  
  out
}

##' @describeIn rprobMat Random conditional distributions
##' @export rcondProbMat
rcondProbMat <- function(n, dim, d, alpha=1, condition) {
  ## if dimension vector shorter than length d, recycle (with warning
  ## if necessary)
  if (missing(d)) d <- length(dim)
  else if (length(dim) < d) dim <- dim * rep.int(1L, d)
  else if (length(dim) > d) stop("More than 'd' dimensions supplied")
  
  if (any(dim < 0)) stop("Dimensions must be non-negative")
  if (any(alpha < 0)) stop("Parameters must be non-negative")
  
  k <- prod(dim)
  out <- rgamma(n*k, alpha, 1)
  dim(out) <- c(n,dim)
  mar <- marginTable(out, c(1,condition+1))
  out <- out/c(mar[patternRepeat0(c(1,condition+1), c(n, dim))])
  dim(out) <- c(n, prod(dim))
  
  class(out) <- "tables"
  attr(out, "tdim") <- dim
  attr(out, "conditional") <- condition
  
  out
}
