#' Calls a G-Wishart sampler
#'
#' This functions takes a matrix X and the output of a previous G-wishart prior sampler and returns a draw from a G-wishart prior in graph form.
#'
#' @param X The data input matrix (D by p)
#' @param bdgraph_obj The output of a previous BDgraph G-wishart sampler
#' @return One draw of a full conditional from a G-wishart prior in list form
#' @export


#find gwish_out
f_gwish_cont = function(X,bdgraph_obj, dfprior){
  sink(tempfile()) 
  on.exit(sink()) 
  out = invisible(force(BDgraph::bdgraph(data = X, g.start = bdgraph_obj, iter = 20, burnin = 0, df.prior = dfprior))) 
  return(out)
}