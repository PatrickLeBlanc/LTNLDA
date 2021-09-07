#' Provides high level summary of model results.
#'
#' This function takes the output of an LTNLDA function as an input, and finds the average subcommunity abundance in decreasing order, as well as the top n ASVs characterizing each subcommunity.  More detailed information in "LTN-LDA" vignette.
#'
#' @param model The output of an LTNLDA model trained on a training set.  
#' @param n We find the n most prevalent ASVs in each subcommunity, with proportions.
#' @return A list with two entries.  Abundance is a vector of the average subcommunity prevalence in all samples in decreasing order. Top_ASVs is a list such tha the kth entry is the the top n ASVs characterizing the kth subcommunity (in decreasing order of average subcommunity prevalence) as well the corresponding ASV prevalencies.
#' @examples
#' data("ps",package = "LTNLDA")
#' K = 2
#' model = LTNLDA(ps,K)
#' sum = Summary(model)
#' @references
#' ADD PAPER INFORMATION ONCE WE KNOW IT
#' @export

Summary = function(model, n = 5){
  
  #check inputs to see if they have been specified or not
  if(is.null(model)){
    cat("Model is not specified. \n")
  } else{
    
    ########
    # Find average abundances and decreasing order
    #recover average ASV-subcommunity distributions
    post_phi_dk = model$Mean_Post_Phi_d
    #recover average abundances
    avg_abund = apply(post_phi_dk,2,mean)
    abund_ord = order(avg_abund,decreasing = TRUE)
    
    abundance = matrix(0,nrow = K, ncol = 1)
    rownames(abundance) = 1:K
    for(k in 1:K){
      abundance[k,1] = avg_abund[abund_ord[k]]
    }
    
    ###########
    # Find top_n ASVs for each subcommunity, in order of average decreasing subcommunity abundnace
    
    #recover average ASV-subcommunity distributions
    post_beta_kv = model$Mean_Post_Beta_k
    #recover phyloseq object
    ps = model$phyloseq
    
    #find the names of the ASVs
    tax = tax_table(ps)
    ASV_names = rownames(tax)
    
    #make a list for the top ASV
    #top_ASV[[k]] contains a matrix listing the top n ASVs in subcommunity k in descending order and with prevlances.
    topn = NULL
    for(k in 1:K){
      new_k = abund_ord[k]
      
      temp = matrix(0,nrow = n, ncol = 1)
      ranked_ASVs = order(post_beta_kv[new_k,],decreasing = TRUE)
      top_ind = ranked_ASVs[1:n]
      temp[1:n,1] = post_beta_kv[new_k,top_ind]
      rownames(temp) = ASV_names[top_ind]
      topn[[k]] = temp
    }
    
    out = list(Abundance = abundance,
               Top_ASVs = topn
    )
    
    return(out)
    
  }
  
}