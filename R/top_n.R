#' Finds the n most prevalent ASVs in each subcommunity for a trained LTN-LDA model.
#'
#' This function takes the output of an LTNLDA function as an input, and finds the top n ASVs characterizing each subcommunity.
#'
#' @param model The output of an LTNLDA model trained on a training set.  
#' @param n We find the n most prevalent ASVs in each subcommunity, with proportions.
#' @return A list such that out[[k]] contains a matrix listing the top n ASVs in subcommunity k in descending order and with the corresponding prevalence.
#' @examples
#' data("ps",package = "LTNLDA")
#' K = 2
#' model = LTNLDA(ps,K)
#' top_n(model)
#' @references
#' ADD PAPER INFORMATION ONCE WE KNOW IT
#' @export


top_n = function(model, n = 5){

  #check inputs to see if they have been specified or not
  if(is.null(model)){
    cat("Model is not specified. \n")
  } else{
    
    #recover average ASV-subcommunity distributions
    post_beta_kv = model$Mean_Post_Beta_k
    #recover phyloseq object
    ps = model$phyloseq
    
    #find the names of the ASVs
    tax = tax_table(ps)
    ASV_names = rownames(tax)
    
    #make a list for the top ASV
    #top_ASV[[k]] contains a matrix listing the top n ASVs in subcommunity k in descending order and with prevlances.
    out = NULL
    for(k in 1:K){
      temp = matrix(0,nrow = n, ncol = 1)
      ranked_ASVs = order(post_beta_kv[k,],decreasing = TRUE)
      top_ind = ranked_ASVs[1:n]
      temp[1:n,1] = post_beta_kv[k,top_ind]
      rownames(temp) = ASV_names[top_ind]
      out[[k]] = temp
    }
    
    return(out)
    
  }
  
}