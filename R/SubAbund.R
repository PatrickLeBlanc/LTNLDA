#' Ranks subcommunities by order of average abundance.
#'
#' This function takes the output of an LTNLDA function as an input, and finds the top n ASVs characterizing each subcommunity.
#'
#' @param model The output of an LTNLDA model trained on a training set.  
#' @return A list such that out[[k]] contains a matrix listing the top n ASVs in subcommunity k in descending order and with the corresponding prevalence.
#' @examples
#' data("ps",package = "LTNLDA")
#' K = 2
#' model = LTNLDA(ps,K)
#' SubAbund(model)
#' @references
#' ADD PAPER INFORMATION ONCE WE KNOW IT
#' @export

SubAbund = function(model){
  
  #check inputs to see if they have been specified or not
  if(is.null(model)){
    cat("Model is not specified. \n")
  } else{
    
    #recover average ASV-subcommunity distributions
    post_phi_dk = model$Mean_Post_Phi_d
    #recover average abundances
    avg_abund = apply(post_phi_dk,2,mean)
    abund_ord = order(avg_abund,decreasing = TRUE)
    
    out = matrix(0,nrow = K, ncol = 1)
    rownames(out) = 1:K
    for(k in 1:K){
      out[k,1] = avg_abund[abund_ord[k]]
    }
    
    return(out)
    
  }
  
}