#' Calls an Inverse-Wishart sampler
#'
#' This functions takes a double v and a scale matrix S to draw from an IW(v,S) distribution.  This is for uses in Cpp code.  
#'
#' @param v The degrees of freedom in an IW(v,S) distribution
#' @param S The scale matrix in an IW(v,S) distribution
#' @return A draw from an IW(v,S) distribution
#' @export


f_iwish = function(v,S){
  return(MCMCpack::riwish(v,S))
}