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
##' @details Returns an object of class \code{tables} consisting of 
##' discrete probability distributions.  Each distribution is assumed to be a 
##' contingency table of dimension \code{dim}, and the probabilities
##' are generated using a Dirichlet distribution with parameters all 
##' equal to \code{alpha}.
##' 
##' @examples
##' dat <- rprobMat(10, c(2,2,2))
##' 
##' @export rprobMat
rprobMat = function(n, dim, d, alpha=1) {
  ## if dimension vector shorter than length d, recycle (with warning
  ## if necessary)
  if (missing(d)) d = length(dim)
  else if (length(dim) < d) dim = dim * rep.int(1L, d)
  else if (length(dim) > d) stop("More than 'd' dimensions supplied")
  
  if (any(dim < 0)) stop("Dimensions must be non-negative")
  if (any(alpha < 0)) stop("Parameters must be non-negative")
  
  k = prod(dim)
  out = matrix(rgamma(n*k, alpha, 1), n, k)
#  out = out/rep(.colSums(out, k, n), each=k)
  out = out/.rowSums(out, n, k)
  
  class(out) = "tables"
  attr(out, "tdim") = dim
#  dim(out) = c(dim, n)
  
  out
}

rcondProbMat <- function(n, dim, d, alpha=1, condition) {
  ## if dimension vector shorter than length d, recycle (with warning
  ## if necessary)
  if (missing(d)) d = length(dim)
  else if (length(dim) < d) dim = dim * rep.int(1L, d)
  else if (length(dim) > d) stop("More than 'd' dimensions supplied")
  
  if (any(dim < 0)) stop("Dimensions must be non-negative")
  if (any(alpha < 0)) stop("Parameters must be non-negative")
  
  k = prod(dim)
  out = rgamma(n*k, alpha, 1)
  dim(out) <- c(n,dim)
  mar <- marginTable(out, c(1,condition+1))
  out = out/c(mar[patternRepeat0(c(1,condition+1), c(n, dim))])
  dim(out) <- c(n, prod(dim))
  
  class(out) = "tables"
  attr(out, "tdim") = dim
  
  out
}

##' Subset object of class tables
##' 
##' Take subset of tables class.
##' @param x object of class \code{tables}
##' @param i indicies of which tables to retain
##' @param j which rows of each table to retain (or if \code{...} not specified, entries )
##' @param ... additional indices up to the dimension of the table
##' @param drop usual logical indicating whether to consolidate margins of the table (doesn't apply to \code{i})
##' @param keep if only one table is specified with \code{i}, should the object output be an object of class \code{tables}?  If not becomes a suitable array.
##' 
##' @details There are two main ways to subset these tables.  In both cases the first index 
##' refers to the tables being selected; one of the methods is to additionally specify all the 
##' indices corresponding to the tables, the other is to only specify a single entry.  
##' For example, \code{x[,1,2,2]} specifies the (1,2,2)th entry of each table; \code{x[,7]} will
##' have the same effect for 2x2x2 tables.
##' 
##' If only one index is specified, then the function behaves just as ordinary subsetting
##' on an array.
##' 
##' @export [.tables
`[.tables` <- function(x, i, j, ..., drop=TRUE, keep=FALSE) {
  mdrop <- missing(drop); mkeep <- missing(keep)
  Nargs <- nargs() - (!mdrop) - (!mkeep)
  if (Nargs <= 2) {
    ## for a single argument (i) just return that entry in the vector
    class(x) <- NULL
    return(x[i])
  }
  else if (Nargs <= 3) {
    ## for two entries treat as a matrix with each table a row
    if (missing(j)) {
      out <- as.matrix(x)[i,,drop=FALSE]
      if (!missing(i) && length(i) == 1 && !keep) {
        dim(out) <- tdim(x)
        dimnames(out) <- tdimnames(x)
      }
      else {
        tdim(out) <- tdim(x)
        tdimnames(out) <- tdimnames(x)
        class(out) <- "tables"
      }
      return(out)
    }
    else {
      out <- as.matrix(x)[i,j,drop=drop]
      return(out)
    }
  }
  else {
    out <- as.array(x)
    out <- out[i, j, ..., drop=drop]
    if (!missing(i) && drop && keep && length(i) <= 1) dim(out) <- c(length(i), dim(out))
    if (missing(i) || length(i) > 1 || keep) {
      class(out) <- "tables"
      tdim(out) <- if (length(dim(out))) dim(out)[-1]
                   else 1L
      dim(out) <- if (length(dim(out))) c(dim(out)[1], length(out)/dim(out)[1])
      else length(out)
    }
    return(out)
  }
  stop()
}

##' Convert tables into array
##' @export as.array.tables
as.array.tables <- function(x) {
   n <- ntables(x)
   dim(x) <- c(n, tdim(x))
#    tdim(x) <- NULL
#     if (!is.null(tdimnames)) {
#      dimnames(x) <- list(seq(n), tdimnames(x))
#       tdimnames(x) <- NULL
#     }
  class(x) <- "array"
  return(x)
}

## maybe need quicker version of as.array()
# regroup <- function(x) {
#   dim(x) <- c(dim(x)[1], tdim(x))
#   class(x) <- "array"
#   return(x)
# }

##' Convert tables into matrix
##' @export as.matrix.tables
as.matrix.tables <- function(x) {
  class(x) <- "matrix"
  attr(x, "tdim") <- NULL
  attr(x, "tdimnames") <- NULL
  return(x)
}

# flatten = function(x) {
#   n = last(dim(x))
#   dim(x) = c(length(x)/n, n)
#   x
# }
# 
# group = function(x) {
#   n = last(dim(x))
#   dim(x) = c(tdim(x), n)
#   x
# }

##' Dimension of distributions over contingency tables
##' @param x an object of class \code{tables}
##' @details The class \code{tables} is used to represent a collection of 
##' multidimentional tables; this function
##' returns the dimension of each table.
##' @return an integer vector of the dimensions
##' @export tdim
##' @export tdim<-
tdim <- function(x) attr(x, "tdim")
`tdim<-` <- function(x, value) {
  attr(x, "tdim") <- value
  x
}
##' Dimension names for distributions over contingency tables
##' @export tdimnames
##' @export tdimnames<-
tdimnames <- function(x) attr(x, "tdimnames")
`tdimnames<-` <- function(x, value) {
  attr(x, "tdimnames") = value
  x
}

cmpfun(`tdim<-`)
cmpfun(tdim)
cmpfun(`tdimnames<-`)
cmpfun(tdimnames)
cmpfun(`[.tables`)
cmpfun(as.matrix.tables)

##' Number of tables
##' @param x an object of class \code{tables}
##' @details Gives the number of tables in an object of class 
##' \code{tables}.
##' @export ntables
ntables <- function(x) {
  dim(x)[1]
}

##' Print tables
##' 
##' Print method for object of class \code{tables}.
##' @param x object of class \code{tables}
##' @param ... arguments to pass to print method for an array
##' @export print.tables
print.tables = function(x, ...) {
  
  if (is.null(dim(x)) || length(dim(x)) == 1) {
    cat("Group of ", length(x), " null numeric tables\n", sep="")
    cat(head(x), sep=", ")
    if (length(x) > 6) cat(", ...\n")
    return(invisible(x))
  }
  n = dim(x)[1L]
  dims = tdim(x)
  k = length(dims)
  
  dim_str = paste(dims, collapse="x", sep="")
  cat("Group of ", n, " numeric tables of dimension ", dim_str, "\n", sep="")
  
  if (n > 0) {cat("First entry:\n")
    y = as.matrix(x)[1,]
    dim(y) = dims
    dimnames(y) = tdimnames(x)
    print.default(y, ...)
  }
  
  invisible(x)
}

# head.tables <- function(x, n=6L, flatten=FALSE, ...) {
#   out = head.matrix(x, n)
#   tdim(out) <- tdim(x)
#   tdimnames(out) <- tdimnames(x)
#   if (is.null(dims)) return(head.default(x, n))
#   N = last(dims)
#   n = min(n, N)
#   y = x[seq(n*prod(dims[-length(dims)]))]
#   if (flatten) attributes(y) <- list(dim=c(length(y)/n,n))
#   else attributes(y) <- list(dim=c(tdim(x),n))
#   y
# }

##' Permute dimensions of tables
##' 
##' Method for permuting indices of tables object.
##' 
##' @param a object of class \code{tables}
##' @param perm permutation of 1,...,k, where each table has k dimensions
##' 
##' @export aperm.tables
aperm.tables <- function(a, perm) {
  mdims <- dim(a)
  dim(a) <- c(ntables(a), tdim(a))
  out <- aperm.default(a, c(1,perm+1))
  dim(out) <- mdims
  
  tdim(out) <- tdim(a)[perm]
  if(!is.null(tdimnames(a))) tdimnames(out) <- tdimnames(a)[perm]
  class(out) <- "tables"
  out
}
