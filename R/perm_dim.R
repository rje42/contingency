##' Permute indices for variable `k`
##' 
##' Currently only works for binary dimensions.
##' 
##' @param x array or related object
##' @param k index to permute
##' @param perm permutation to perform
##' @param ... other arguments (not currently used)
##' 
##' @details Permutes the levels of one variable 
##' according to the permutation given in `perm`.
##' Can be applied to matrices, arrays or tables.
##' 
##' @return A permuted `array` or `tables` object.
##' 
##' @export 
perm_dim <- function(x, k, perm, ...) {
  UseMethod("perm_dim")
}

##' @export
perm_dim.default <- function(x, k, perm, ...) {
  y <- x
  
  ## check this is a valid permutation, or just shift one along
  if (missing(perm)) {
    perm <- c(seq_len(dim(x)[k])[-1], 1)
  }
  else if (!isTRUE(all.equal(sort.int(perm), seq_len(dim(x)[k])))) {
    stop("Not a valid permutation")
  }
  
  ## implement the permutation
  for (i in seq_len(dim(x)[k])) {
    subtable(y,k,i) <- subtable(x,k,perm[i])
  }

  return(y)
}

##' @export
perm_dim.tables <- function(x, k, perm, ...) {
  
  if (missing(perm)) {
    perm <- c(seq_len(tdim(x)[k])[-1], 1)
  }

  ## call default method to perform permutation  
  out <- perm_dim.default(x=as.array.tables(x), k=k+1, perm=perm)
  out <- as_tables(out, tdim=tdim(x), rev=is_rev(x))
  # dim(out) <- dim(x)
  # tdim(out) <- tdim(x)
  
  return(out)
}
