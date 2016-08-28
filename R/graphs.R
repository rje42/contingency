##' Generate probability distributions according to a graphical model
##' 
##' Generate tables which factorize according to a directed acyclic graph.
##' 
##' @param graph object of class \code{mixedgraph} giving a directed acyclic graph
##' @param n number of distributions to generate
##' @param alpha parameter controlling Dirichlet distributions
##' @param dim integer vector of dimensions for the distribution
##' 
##' @export rDAGmodel
##' 
##' @examples 
##' library(MixedGraphs)
##' gr <- mixedgraph(5, edges=list(directed=list(c(1,2),c(1,3),c(2,4),c(3,5),c(4,5))))
##' rDAGmodel(gr, 1)
rDAGmodel <- function(graph, n, dim=2L, alpha=1) {
  require(MixedGraphs)
  ## if not currently in mixedgraph format, try coercing it to be so
  if (!is.mixedgraph(graph)) graph = convert(graph, format="mixedgraph")
  d <- length(graph$v)

  ## if dimension vector shorter than length d, recycle 
  ## (with warning if necessary)
  if (length(dim) < d) dim = dim * rep.int(1L, d)
  else if (length(dim) > d) stop("More than 'd' dimensions supplied")
  
  ord <- topologicalOrder(graph)  
  out <- matrix(1, nrow=n, ncol=prod(dim))
  
  for (v in seq_along(ord)) {
    pav <- pa(graph, v)
    tmp <- rcondProbMat(n, dim = dim[c(v,pav)], alpha=alpha, condition=seq_along(pav)+1)
    tmp <- c(tmp[,patternRepeat0(c(v,pav),dim,keep.order=TRUE)])
    out <- out*c(tmp)
  }
  class(out) <- "tables"
  tdim(out) <- dim

  out
}


