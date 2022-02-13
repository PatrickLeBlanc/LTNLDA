#' Calls a G-Wishart sampler
#'
#' This functions takes a matrix X and returns a draw from a G-wishart prior in graph form.
#'
#' @param X The data input matrix (D by p)
#' @return One draw of a full conditional from a G-wishart prior in list form
#' @export


#find gwish_out
f_gwish = function(X,dfprior){
  sink(tempfile()) 
  on.exit(sink()) 
  out = invisible(force(BDgraph::bdgraph(data = X,  iter = 20, burnin = 0, df.prior = dfprior))) 
  return(out)
}