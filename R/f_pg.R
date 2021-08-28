#' Calls a Polya-Gamma Sampler
#'
#' This functions takes two numbers b,c as inputs and returns a draw from a PG(b,c) distribution and is written to allow calling of the pgdraw function in the Rcpp code.  For usage in R code, it's better to use the pgdraw package directly.
#'
#' @param b The first parameter in a PG(b,c) distribution
#' @param c The second parameter in a PG(b,c) distribution
#' @return A draw from a PG(b,c) distribution
#' @references 
#' Nicholas G. Polson, James G. Scott, and Jesse Windle.  Bayesian Inference for logistic models using polya-gamma latent variables.  Journal of the American Statistical Association.  108(504):1339-1349. 2013.
#' @export


f_pg = function(b,c){
  return(pgdraw::pgdraw(b,c))
}