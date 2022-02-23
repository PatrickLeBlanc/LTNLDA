#' Calls a G-Wishart sampler
#'
#' This functions takes an adjaceny matrix G, a degree of freedom m, and a scale matrix M, and returns a draw form a W_G(m,M) distribution
#'
#' @param G An upper-triangular adjaceny matrix
#' @param m The degree of freedom in  a G-Wishart distribution
#' @param M The scale matrix in a G-Wishart distribution
#' @return One draw from a W_G(m,M) distribution
#' @export

#Draw from an exact g-wishart sampler
f_gwish = function(G,m,M){
  return(BDgraph::rgwish(n = 1, adj = G, b = m, D = M))
}