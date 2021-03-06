% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f_pg.R
\name{f_pg}
\alias{f_pg}
\title{Calls a Polya-Gamma Sampler}
\usage{
f_pg(b, c)
}
\arguments{
\item{b}{The first parameter in a PG(b,c) distribution}

\item{c}{The second parameter in a PG(b,c) distribution}
}
\value{
A draw from a PG(b,c) distribution
}
\description{
This functions takes two numbers b,c as inputs and returns a draw from a PG(b,c) distribution and is written to allow calling of the pgdraw function in the Rcpp code.  For usage in R code, it's better to use the pgdraw package directly.
}
\references{
Nicholas G. Polson, James G. Scott, and Jesse Windle.  Bayesian Inference for logistic models using polya-gamma latent variables.  Journal of the American Statistical Association.  108(504):1339-1349. 2013.
}
