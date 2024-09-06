##' Test if tables object is reversed
##' 
##' Checks if the rows of the underlying matrix represent each table or each entry in a table
##' 
##' @param x object to be tested
##' 
##' @export
is_rev <- function (x) {
  if (!is_tables(x)) return(FALSE)
  
  if (is.null(attr(x, "rev"))) return(FALSE)
  
  return(attr(x, "rev"))
}

