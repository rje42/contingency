##' Wrappers for Dirichlet distribution over a table
##' 
##' @param x `tables` object of observations
##' @param n number of samples
##' @param alpha table containing parameters
##' @param rev logical: should output move through each table fastest?
##' 
##' @details
##' This function obtains the Dirichlet density over a contingency table 
##' structure.  In other words, suppose that we have a matrix observation 
##' \eqn{x = (x_{ij})} where \eqn{\sum_{i,j} x_{i,j} = 1}.  Then we might choose
##' to model the vector \eqn{x} as having a Dirichlet distribution, with weights
##' \eqn{\alpha_{ij}}.
##' 
##' If `alpha` is a scalar in `dtab_dir` then it is applied to every entry in `x`.
##' 
##' @name tab_dir
NULL

##' @describeIn tab_dir density function
##' @export
dtab_dir <- function (x, alpha, log = FALSE) {
  contingency::capply(x, function (z) ddirichlet(c(z), alpha=c(alpha), log=log))
}

##' @describeIn tab_dir sampling function
##' @export
rtab_dir <- function (n, alpha, rev=FALSE) {
  x <- rdirichlet(n, c(alpha))
  x <- as_tables(x, tdim=dim(alpha), rev=rev)
  
  return(x)
}
