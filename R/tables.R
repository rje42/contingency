##' Subset object of class tables
##' 
##' Take subset of tables class.
##' @param x object of class `tables`
##' @param i indices of which tables to retain
##' @param j which rows of each table to retain (or if `...` not specified, entries )
##' @param ... additional indices up to the dimension of the table
##' @param drop usual logical indicating whether to consolidate margins of the table (doesn't apply to `i`)
##' @param keep if only one table is specified with `i`, should the object output be an object of class `tables`?  If not becomes a suitable array.
##' 
##' @details There are two main ways to subset these tables.  In both cases the first index 
##' refers to the tables being selected; one of the methods is to additionally specify all the 
##' indices corresponding to the tables, the other is to only specify a single entry.  
##' For example, `x[,1,2,2]` specifies the (1,2,2)th entry of each table; `x[,7]` will
##' have the same effect for 2x2x2 tables.
##' 
##' If only one index is specified, then the function behaves just as ordinary 
##' subsetting on an array, though a reversed table is modified to have its 
##' first dimension index the tables.
##' 
##' @return A tables object over the specific entries and values selected.
##' 
##' @examples 
##' x <- rprobMat(n=10, rep(2,3))
##' x[1,]
##' x[,1,1:2,1]
##' x[,1,1:2,1,drop=FALSE]
##' 
##' @export
`[.tables` <- function(x, i, j, ..., drop=TRUE, keep=FALSE) {
  mdrop <- missing(drop); mkeep <- missing(keep)
  Nargs <- nargs() - (!mdrop) - (!mkeep)

  rev <- is_rev(x)
  if (rev) x <- reverse(x)
    
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
        attr(out, "conditional") <- attr(x, "conditional")
        if (rev) out <- reverse(out)
      }
      return(out)
    }
    else {
      if (rev) out <- as.matrix(x)[j,i,drop=drop]
      else out <- as.matrix(x)[i,j,drop=drop]
      return(out)
    }
  }
  else {
    out <- as.array(x)
    attr(out, "conditional") <- attr(x, "conditional")
    out <- out[i, j, ..., drop=drop]
    if (!missing(i) && drop && keep && length(i) <= 1) dim(out) <- c(length(i), dim(out))
    if (missing(i) || length(i) > 1 || keep) {
      class(out) <- "tables"
      tdim(out) <- if (length(dim(out))) dim(out)[-1]
      else 1L
      dim(out) <- if (length(dim(out))) c(dim(out)[1], length(out)/dim(out)[1])
      else length(out)
    }
    if (rev) {
      if (is_tables(out)) out <- reverse(out)
      else if (!is.null(dim(out))) {
        nd <- length(dim(out))
        if (nd > 1) out <- aperm(out, c(seq_len(nd-1)+1, 1))
      }
    }
    return(out)
  }
  stop()
}

##' Convert tables into array
##' 
##' @param x `tables` object
##' @param ... other arguments
##' 
##' @return An `array` object
##' 
##' @method as.array tables
##' @export 
as.array.tables <- function(x, ...) {
  n <- ntables(x)
  dim(x) <- c(n, tdim(x))
  #    tdim(x) <- NULL
  #     if (!is.null(tdimnames)) {
  #      dimnames(x) <- list(seq(n), tdimnames(x))
  #       tdimnames(x) <- NULL
  #     }
  class(x) <- "array"
  wh <- match(c("tdim", "tdimnames", "conditional"), names(attributes(x)), nomatch = 0L)
  if (any(wh) > 0) attributes(x) <- attributes(x)[-wh]
  
  return(x)
}

## maybe need quicker version of as.array()
# regroup <- function(x) {
#   dim(x) <- c(dim(x)[1], tdim(x))
#   class(x) <- "array"
#   return(x)
# }

##' Convert tables into matrix
##' 
##' @param x `tables` object
##' @param ... other arguments
##' 
##' @return A `matrix` object
##' 
##' @method as.matrix tables
##' @export 
as.matrix.tables <- function(x, ...) {
  class(x) <- "matrix"
  wh <- match(c("tdim", "tdimnames", "conditional"), names(attributes(x)), nomatch = 0L)
  if (any(wh) > 0) attributes(x) <- attributes(x)[-wh]
  return(x)
}

# flatten <- function(x) {
#   n <- last(dim(x))
#   dim(x) <- c(length(x)/n, n)
#   x
# }
# 
# group <- function(x) {
#   n <- last(dim(x))
#   dim(x) <- c(tdim(x), n)
#   x
# }

##' Dimension of distributions over contingency tables
##' 
##' @param x an object of class `tables`
##' 
##' @details The class `tables` is used to represent a collection of 
##' multidimentional tables; this function
##' returns the dimension of each table.
##' 
##' @return an integer vector of the dimensions
##' 
##' @export tdim
tdim <- function(x) attr(x, "tdim")
##' @describeIn tdim assign tables dimension
##' @param value value to set parameters to
##' @return the `tables` object inputted with the new dimensions
##' @export tdim<-
`tdim<-` <- function(x, value) {
  attr(x, "tdim") <- value
  x
}

##' Dimension names for distributions over contingency tables
##' 
##' @param x `tables` object
##' 
##' @export tdimnames
tdimnames <- function(x) attr(x, "tdimnames")
##' @describeIn tdimnames assign dimension names
##' @param value value to set dimension names to
##' @return the `tables` object inputted with the new dimension names
##' @export tdimnames<-
`tdimnames<-` <- function(x, value) {
  attr(x, "tdimnames") <- value
  x
}

# cmpfun(`tdim<-`)
# cmpfun(tdim)
# cmpfun(`tdimnames<-`)
# cmpfun(tdimnames)
# cmpfun(`[.tables`)
# cmpfun(as.matrix.tables)

##' Number of tables
##' @param x an object of class `tables`
##' @details Gives the number of tables in an object of class 
##' `tables`.
##' @return An integer.
##' @export ntables
ntables <- function(x) {
  if (is_rev(x)) dim(x)[2]
  else dim(x)[1]
}

##' Print tables
##' 
##' Print method for object of class `tables`.
##' @param x object of class `tables`
##' @param ... arguments to pass to print method for an array
##' 
##' @return The input provided (invisibly).
##' 
##' @method print tables
##' @export 
print.tables <- function(x, ...) {
  
  if (is.null(dim(x)) || length(dim(x)) == 1) {
    cat("Group of ", length(x), " null numeric tables\n", sep="")
    cat(head(x), sep=", ")
    if (length(x) > 6) cat(", ...\n")
    return(invisible(x))
  }
  
  rev <- is_rev(x)
  
  if (rev) n <- dim(x)[2L] 
  else n <- dim(x)[1L]
  dims <- tdim(x)
  k <- length(dims)
  
  dim_str <- paste(dims, collapse="x", sep="")
  cat("Group of ", n, ifelse(rev, " reversed",""), " numeric tables of dimension ", dim_str, "\n", sep="")
  
  if (n > 0) {cat("First entry:\n")
    if (rev) y <- x[,1L]
    else y <- x[1L,]
    dim(y) <- dims
    dimnames(y) <- tdimnames(x)
    print.default(y, ...)
  }
  
  invisible(x)
}

# head.tables <- function(x, n=6L, flatten=FALSE, ...) {
#   out <- head.matrix(x, n)
#   tdim(out) <- tdim(x)
#   tdimnames(out) <- tdimnames(x)
#   if (is.null(dims)) return(head.default(x, n))
#   N <- last(dims)
#   n <- min(n, N)
#   y <- x[seq(n*prod(dims[-length(dims)]))]
#   if (flatten) attributes(y) <- list(dim=c(length(y)/n,n))
#   else attributes(y) <- list(dim=c(tdim(x),n))
#   y
# }

##' Permute dimensions of tables
##' 
##' Method for permuting indices of tables object.
##' 
##' @param a object of class `tables`
##' @param perm permutation of 1,...,k, where each table has k dimensions
##' @param ... other arguments to methods
##' 
##' @method aperm tables
##' 
##' @return A permuted `tables` object.
##' 
##' @export
aperm.tables <- function(a, perm, ...) {
  mdims <- dim(a)
  if (is_rev(a)) {
    dim(a) <- c(tdim(a), ntables(a))
    out <- aperm.default(a, c(perm,length(dim(a))))
  }
  else {
    dim(a) <- c(ntables(a), tdim(a))
    out <- aperm.default(a, c(1,perm+1))
  }
  dim(out) <- mdims
  
  attr(out, "conditional") <- order(perm)[attr(out, "conditional")]
  tdim(out) <- tdim(a)[perm]
  if(!is.null(tdimnames(a))) tdimnames(out) <- tdimnames(a)[perm]
  class(out) <- "tables"
  out
}

##' Test object is a collection of tables
##' 
##' Checks if `"tables"` is in the list of class attributes.
##' 
##' @param x object to be tested
##' 
##' 
##' @export
is_tables <- function (x) {
  return("tables" %in% class(x))
}

##' As tables
##' 
##' @param x array or matrix object
##' @param tdim dimensions for each table
##' @param conditional integer vector of indices that are conditional
##' @param rev logical: should output move through each table fastest?
##' @param ... other arguments for methods
##' 
##' @details Transforms a vector, matrix or array into a tables object with 
##' the specified dimensions.  Note that if `rev = TRUE` then the final dimension
##' is (by default) used to index the tables objects, and otherwise the first 
##' dimension..
##' 
##' @return A `tables` object.
##' 
##' @export
as_tables <- function(x, tdim, conditional, rev=FALSE, ...) {
  UseMethod("as_tables")
}

##' @method as_tables default
##' @export
as_tables.default <- function(x, tdim, conditional, rev=FALSE, ...) {
  class(x) <- "tables"
  if (missing(conditional)) attr(x, "conditional") <- integer(0)
  else attr(x, "conditional") <- conditional
  
  if (missing(tdim)) {
    contingency::tdim(x) <- length(x)
    if (rev) dim(x) <- c(length(x), 1)
    else dim(x) <- c(1, length(x))
  }
  else {
    contingency::tdim(x) <- tdim
    N <- length(x)/prod(tdim)
    if (N != ceiling(N)) warning("Supplied dimensions don't give a whole number of tables")
    if (rev) dim(x) <- c(prod(tdim), ceiling(N))
    else dim(x) <- c(ceiling(N), prod(tdim))
  }
  
  return(x)
}

##' @method as_tables matrix
##' @export
as_tables.matrix <- function(x, tdim, conditional, rev=FALSE,...) {
  class(x) <- "tables"
  if (missing(conditional)) attr(x, "conditional") <- integer(0)
  else attr(x, "conditional") <- conditional
  
  if (missing(tdim)) {
    contingency::tdim(x) <- dim(x)
    if (rev) dim(x) <- c(prod(dim(x)), 1)
    else dim(x) <- c(1, prod(dim(x)))
  }
  else {
    contingency::tdim(x) <- tdim
    N <- length(x)/prod(tdim)
    if (N != ceiling(N)) warning("Supplied dimensions don't give a whole number of tables")
    if (rev) dim(x) <- c(prod(tdim), ceiling(N))
    else dim(x) <- c(ceiling(N), prod(tdim))
  }
  
  return(x)
}

##' @method as_tables array
##' @export
as_tables.array <- function(x, tdim, conditional, rev=FALSE, ...) {
  class(x) <- "tables"
  if (missing(conditional)) attr(x, "conditional") <- integer(0)
  else attr(x, "conditional") <- conditional
  
  if (missing(tdim)) {
    contingency::tdim(x) <- dim(x)
    if (rev) dim(x) <- c(prod(dim(x)), 1)
    else dim(x) <- c(1, prod(dim(x)))
  }
  else {
    contingency::tdim(x) <- tdim
    N <- length(x)/prod(tdim)
    if (N != ceiling(N)) warning("Supplied dimensions don't give a whole number of tables")
    if (rev) dim(x) <- c(prod(tdim), ceiling(N))
    else dim(x) <- c(ceiling(N), prod(tdim))
  }
  
  return(x)
}

##' Create blank tables
##' 
##' @param n number of tables
##' @param tdim dimension of each table
##' @param rev logical: should output move through each table fastest?
##' 
##' @export
tables <- function(n, tdim, rev=FALSE) {
  as_tables.array(array(1, dim=c(n, tdim)), tdim=tdim, rev=rev)
}

##' Bind tables of the same dimension
##' 
##' @param x a `tables object`
##' @param ... further `tables` objects with the same `tdim` attributes
##' @param rev logical: should output move through each table fastest?
##' 
##' @export
tbind <- function (x, ..., rev=FALSE) {
  lst <- list(x, ...)
  
  crev <- sapply(lst, is_rev)
  lst[crev] <- lapply(lst[crev], reverse)
  
  dims <- sapply(lst, tdim)
  if (is.list(dims)) stop("Different numbers of dimensions")
  else if (is.matrix(dims)) {
    if (any(apply(dims, 1, sd) > 0)) stop("Different dimensions")
  }
  else if (is.vector(dims)) {
    if (sd(dims) > 0) stop("Different dimensions")
  }
  else stop("Object not expected")
  
  ## now combine objects and return
  out <- do.call(rbind, lst)
  out <- as_tables(out, tdim=tdim(x))
  
  if (rev) out <- reverse(out)
  
  return(out)
}

reverse <- function (x) {
  if (!is_tables(x)) stop("'reverse' only applies to 'tables' objects")
  
  x <- t(x)
  if (is.null(attr(x, "rev"))) attr(x, "rev") <- TRUE
  else attr(x, "rev") <- !attr(x, "rev")
  
  return(x)
}

# as.tables.data.frame <-
