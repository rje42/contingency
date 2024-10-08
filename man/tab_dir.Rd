% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tab_dir.R
\name{tab_dir}
\alias{tab_dir}
\alias{dtab_dir}
\alias{rtab_dir}
\title{Wrappers for Dirichlet distribution over a table}
\usage{
dtab_dir(x, alpha, log = FALSE)

rtab_dir(n, alpha, rev = FALSE)
}
\arguments{
\item{x}{\code{tables} object of observations}

\item{alpha}{table containing parameters}

\item{n}{number of samples}

\item{rev}{logical: should output move through each table fastest?}
}
\description{
Wrappers for Dirichlet distribution over a table
}
\details{
This function obtains the Dirichlet density over a contingency table
structure.  In other words, suppose that we have a matrix observation
\eqn{x = (x_{ij})} where \eqn{\sum_{i,j} x_{i,j} = 1}.  Then we might choose
to model the vector \eqn{x} as having a Dirichlet distribution, with weights
\eqn{\alpha_{ij}}.

If \code{alpha} is a scalar in \code{dtab_dir} then it is applied to every entry in \code{x}.
}
\section{Functions}{
\itemize{
\item \code{dtab_dir()}: density function

\item \code{rtab_dir()}: sampling function

}}
