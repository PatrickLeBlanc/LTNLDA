#' Calls a Bayesian GLasso sampler
#'
#' This functions takes a matrix X and two numbers r and q as inputs and returns a draw from the full conditional of a Bayesian Graphical LASSO distribution.
#'
#' @param X The input matrix
#' @param r The shape hyperparameter in a gamma distribution
#' @param q The rate hyperparameter in a gamma distribution
#' @return One iterate of a Bayesian GLasso full conditional
#' @references 
#' Hao Wang.  Bayesian Graphical Lasso Models and Efficient Posterior Computation.  Bayesian Anal. 7(4):867-886.  December 2012
#' @export


f_bglasso = function(X,r,q){
  return(BayesianGLasso::blockGLasso(X,iterations = 1, burnIn = 0,lambdaPriora = r, lambdaPriorb = q)) 
}
