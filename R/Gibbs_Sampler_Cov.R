#' Runs an LTN-LDA Gibbs sampler.  Experimental, trying out a new covariance structure.  Run at own risk.  
#'
#' This function takes a phyloseq object and the modelled number of subcommunities as inputs. It runs a collapsed blocked Gibbs sampler for LTNLDA, and returns a list containing posterior mean estimates for some parameters, Markov Chains for all meaningful parameters, and the phyloseq object.  For a detailed example see the vignette "LTN-LDA".
#'
#' @param ps A phyloseq object containing an otu_table() and a phy_tree() with an edge matrix.  That is, otu_table(ps) and phy_tree(ps)$edge both exist.  Further, the otu_table(ps) should have taxa corresponding to rows adn samples corresponding to columns.
#' @param K An integer specifying the number of modeled subcommunities.
#' @param C An integer specifying the threshold controlling cross-sample heterogeneity.  The default value is 5.  Using the default value may result in unreliable inference.
#' @param iterations The number of iterations to record.  Default value is 1000.
#' @param burnin The number of burnin iterations to run before recording values.  The default value is 10000.
#' @param thin The amount by which we thin the chain.  A value of X means that 1 every X values is recorded after the burnin.  The default value is 10.
#' @param alpha A double specifying the prior on the subcommunity-sample proportions.  The default value is 1.
#' @param a_U A double controlling the scale parameter in the inverse-Gamma distribution for the upper covariance values.  The default value is 10^4.
#' @param b_U A double controlling the rate parameter in the inverse-Gamma distribution for the upper covariance values.  The default value is 10.
#' @param a_L A double such that the degrees of freedom for the G-Wishart prior for the lower precision matrix is a_L(p_L+2).  The default value is 5.
#' @param b_L A double such that the scale matrix for the G-Wishart prior for the lower precision matrix is b_L(p+2)diag(p_L.  The default value is 5.
#' @param g_prior A double specifying the prior probability that there is an edge between two nodes in the lower covariance matrix.  The default value is 1/4.
#' @param Lambda A matrix specifying a covariance prior for the mu_k.  The default value is diag(V) where V is the number of leaves.
#' @return 
#' A list with 8 entries.  
#' Mean_Post_Phi_d contains the posterior mean estimate for the subcommunity-sample distributions phi_d. 
#' Mean_Post_Beta_kd contains the posterior mean estimate for the sample and subcommunity specific multinomial distributions beta_{k,d}.
#' Mean_Post_Beta_k contains the posterior mean estimates for the average subcommunity multinomial distributions beta_k. 
#' Mean_Post_G_L contains the posterior mean of the edge matrix denoting the probability edges are conditionally dependent for each subcomunity.
#' Chain_Phi contains the Markov Chains for the phi_d. 
#' Chain_Psi contains the Markov Chains for the psi_{p,d,k}. 
#' Chain_Mu contains the Markov Chains for the mu_k. 
#' Chain_Sigma contains the Markov Chains for the Sigma_k. 
#' phyloseq contains the phyloseq object the Gibbs sampler ran on.  
#' @examples
#' data("ps",package = "LTNLDA")
#' K = 2
#' model = LTNLDA(ps,K)
#' @references
#' ADD PAPER INFORMATION ONCE WE KNOW IT
#' @export


LTNLDA_cov = function(ps, K, C = 5,
                  iterations = 1000, burnin = 10000, thin = 10,
                  alpha = 1, a_L = 5, b_L = 5, a_U = 10^4, b_U = 10, 
                  g_prior = 1/4,
                  Lambda = NULL){
  
  #check inputs to see if they have been specified or not
  if(is.null(ps)){
    cat("phyloseq object is not specified. \n")
  } else{
    if(is.null(K)){
      cat("K is not specified. \n")
    } else{
      
      cat("Processing data. \n")
      
      #extract structures from phyloseq
      tree.edge = phyloseq::phy_tree(ps)$edge
      dtm = phyloseq::otu_table(ps)
      
      ####################################################
      # Convert Tree Edge Matrix into useable structures #
      ####################################################
      
      #find characteristics of tree
      #find the root node
      root = setdiff(tree.edge[,1],tree.edge[,2]) 
      
      #find the internal nodes and reorder
      internal_nodes = unique(tree.edge[,1])
      internal_nodes = sort(internal_nodes)
      internal_nodes_C = internal_nodes - 1
      
      #find the maximum of the internal nodes
      A = max(internal_nodes)
      
      #find the set of leaves
      leaves = setdiff(tree.edge[,2],tree.edge[,1])
      
      #find the number of leaves
      V = length(leaves)
      
      #Generate some tree data structures
      
      
      #descendants[[i]] contains the two immediate descendants of node i
      descendants = NULL
      descendants_mat = matrix(0,ncol=2,nrow=max(tree.edge))
      for(i in 1:max(tree.edge)){
        if (sum(which(tree.edge[,1]==i))>0){
          descendants[[i]] = tree.edge[which(tree.edge[,1] == i),2]
          descendants_mat[i,] = descendants[[i]]
        }
      }
      descendants_mat_C = descendants_mat - 1
      
      
      #parents[i] contains the parent of node i
      parents = NULL
      for (i in 1:max(tree.edge)){
        if (sum(which(tree.edge[,2]==i))>0){
          parents[i] = tree.edge[which(tree.edge[,2]==i),1]
        }
      }
      
      #ancestors[[i]] contains all of the ancestors of node i
      ancestors = NULL
      ancestors_C = NULL
      for (i in 1:max(tree.edge)){
        up=NULL
        parent = parents[i]
        while (is.na(parents[parent])==FALSE){
          up = c(up,parent)
          parent = parents[parent]
        }
        ancestors[[i]] = c(up,parent) #Adds the root of the tree as well
        ancestors_C[[i]] = ancestors[[i]]-1
      }
      
      #layers[[i]] containts the nodes in the i^th layer of the tree
      #layers[[1]] is the root node, layers[[2]] is the root nodes children, etc
      layers = NULL
      #initialize layer 2 and 1
      #have to do an n-tuple (n>1) first, o/w this explodes
      layers[[2]] = descendants[[root]]
      layers[[1]] = root
      for (i in 3:max(tree.edge)){
        descend = NULL
        for (j in 1:length(layers[[i-1]])){
          descend = c(descend,descendants[[layers[[i-1]][j]]])
        }
        if ((sum(descend)>0)==TRUE){
          layers[[i]] = descend
        } else{
          break
        }
      }
      
      #left_leaves[[node]] contains the leaves left-descended from node
      #right_leaves[[node]] contains the leaves right-descended from node
      left_leaves = NULL
      right_leaves = NULL
      left_leaves[[max(internal_nodes)+1]] = rep(0,5)
      right_leaves[[max(internal_nodes)+1]] = rep(0,5)
      for (node in internal_nodes){
        left_descend = NULL
        right_descend = NULL
        
        descend = descendants[[node]]
        left = descend[1]
        right = descend[2]
        
        #if the descendant is a leaf we can termiante
        if((left %in% leaves)==TRUE){
          left_descend = left
        } else {
          #cycle through all of the leaves and see which are left descendants
          for (nodes in leaves){
            if((left %in% ancestors[[nodes]])==TRUE){
              left_descend = c(left_descend,nodes)
            }
          }
        }
        left_leaves[[node]] = left_descend
        
        #if the descendant is a leaf we can termiante
        if((right %in% leaves)==TRUE){
          right_descend = right
        } else {
          #cycle through all of the leaves and see which are right descendants
          for (nodes in leaves){
            if((right %in% ancestors[[nodes]])==TRUE){
              right_descend = c(right_descend,nodes)
            }
          }
        }
        right_leaves[[node]] = right_descend
      }
      left_leaves[[max(internal_nodes)+1]] = NULL
      right_leaves[[max(internal_nodes)+1]] = NULL
      
      
      
      #need to find, for each leaf
      # the nodes which have to succeed
      # the nodes which have to fail
      #in order for the leaf to be selected
      
      #leaf_success[[leaf]] contains the nodes from which leaf is left-descended
      leaf_success = NULL
      leaf_success_C = NULL
      leaf_success[[max(leaves)+1]] = rep(0,5)
      leaf_success_C[[max(leaves)+1]] = rep(0,5)
      for (leaf in leaves){
        node_list=NULL
        node_list = c(leaf,ancestors[[leaf]])
        successes = NULL
        for (node in ancestors[[leaf]]){
          if ((descendants[[node]][1] %in% node_list)==TRUE){
            successes = c(successes,node)
          }
        }
        if (is.null(successes)==FALSE){
          leaf_success[[leaf]] = successes
          leaf_success_C[[leaf]] = successes - 1
        } else {
          leaf_success_C[[leaf]] = -1
        }
      }
      leaf_success[[max(leaves)+1]]  = NULL
      leaf_success_C[[max(leaves)+1]]  = NULL
      
      #leaf_failures[[leaf]] contains the nodes from which leaf is right-descended
      leaf_failures = NULL
      leaf_failures[[max(leaves)+1]] = rep(0,5)
      leaf_failures_C = NULL
      leaf_failures_C[[max(leaves)+1]] = rep(0,5)
      for (leaf in leaves){
        node_list=NULL
        node_list = c(leaf,ancestors[[leaf]])
        failures = NULL
        for (node in ancestors[[leaf]]){
          if ((descendants[[node]][2] %in% node_list)==TRUE){
            failures = c(failures,node)
          }
        }
        if (is.null(failures)==FALSE){
          leaf_failures[[leaf]] = failures
          leaf_failures_C[[leaf]] = failures - 1
        } else {
          leaf_failures_C[[leaf]] = -1
        }
      }
      leaf_failures[[max(leaves)+1]] = NULL
      leaf_failures_C[[max(leaves)+1]] = NULL
      
      #come up with a mapping from internal nodes to 1:p
      #node_map[internal_nodes] in {1,2,dots,p}
      p = length(internal_nodes)
      node_map = rep(0,A)
      for(x in 1:p){
        node_map[internal_nodes[x]] = x
      }
      
      
      #find the number of leaves descended from each node
      #num_leaves[a] is the number of leaves descended from node a
      num_leaves = rep(0,p)
      for(x in 1:p){
        num_leaves[x] = length(c(left_leaves[[internal_nodes[x]]], right_leaves[[internal_nodes[x]]]))
      }
      
      #make data structure recording which nodes belong in the upper part of the tree and which belong in the lower part
      #U_nodes list the nodes in the upper matrix; p_U is the number of such nodes
      #L_nodes list the nodes in the upper matrix; p_L is the number of such nodes
      U_nodes = which(num_leaves >= C)
      U_nodes_C = U_nodes - 1
      p_U = length(U_nodes)
      L_nodes = which(num_leaves < C)
      L_nodes_C = L_nodes - 1
      p_L = length(L_nodes)
      
      #define scale matrices for Inverse-Wishart covariance priors for upper and lower matrices
      #initialized to diagonal matrices
      Phi_U = diag(p_U)
      Phi_L = diag(p_L)
      
      #scale hyperprior for b_L
      b_L = (p_L+2)*b_L
      
      #check if Lambda has been specified --- if not, specify as identity matrix
      if(is.null(Lambda)){
        Lambda = diag(p)
      }
      
      
      #######################################################
      # Convert dtm count data to list of vectors of tokens #
      #######################################################
      
      #find number of samples
      D = ncol(dtm) 
      
      #convert dtm count data to list of vectors of tokens
      docs = NULL
      docs[[D+1]] = rep(0,10000)
      docs_C = NULL
      for (d in 1:D){ #for each document
        doc = NULL
        
        for(v in 1:V){
          doc = c(doc,rep(v,dtm[v,d]))
        }
        
        doc = doc[sample(1:length(doc),replace = FALSE)]
        docs[[d]] = doc
        docs_C[[d]] = docs[[d]] - 1
      }
      docs[[D+1]] = NULL
      
      ################################################
      # Initialize Data structures for Gibbs Sampler #
      ################################################
      
      # Initialize token-subcommunity assignments and counts
      
      #1 Initialize subocmmunity assignments to sequencing in each sample 
      ta = lapply(docs, function(x) rep(0, length(x))) # initialize topic assignment list
      ta_C = lapply(docs, function(x) rep(0, length(x))) 
      #2 Generate word-topic (ASV-subcommunity) count matrix.
      wt = matrix(0, K, V) # wt[k,v] is the count of word v assigned to topic k
      #3 Initialize document-topic (sample-subcommunity) count matrix
      dt = matrix(0, length(docs), K) #dt[d,k] is count of topic k in document d
      #4 Initialize counts of words-by-document-by topic
      # (ASVs-by-sample-by-subcommunity)
      wc_dwt = array(0,dim=c(D,V,K)) 
      #wc_dwt[d,w,k] is count of word w assigned to topic k in document d
      #5 Initialie counts at each node by document and topic (sample and subcommunity)
      nc_dnt =  array(0,dim=c(D,max(tree.edge),K)) #node-count-by-document-by-topic
      #nc_dnt[d,n,t] is count of words in document d assigned to topic k descended from node n
      
      for(d in 1:length(docs)){ # for each document
        
        for(w in 1:length(docs[[d]])){ # for each token in document d
          ta[[d]][w] = sample(1:K, 1) # randomly assign topic to token
          ta_C[[d]][w] = ta[[d]][w]-1
          ti = ta[[d]][w] # topic index
          wi = docs[[d]][w] # word associated with token w
          wt[ti,wi] = wt[ti,wi]+1 # update word-topic count matrix     
        }
        
        #Use subcommunity-assignment vector to:
        #   find sample-subcommunity (dt) counts
        #   find ASVs-by-sample-by-subcommunity (wc_dwt) counts
        for(k in 1:K){ # for each topic k
          dt[d,k] = sum(ta[[d]]==k) # count tokens in document d assigned to topic t 
          
          ta.temp = docs[[d]][which(ta[[d]]==k)]
          for(w in 1:length(ta.temp)){ # for each word
            ti = k # topic index
            wi = ta.temp[w] # wordID for token w
            wc_dwt[d,wi,ti] = wc_dwt[d,wi,ti]+1 # update word-topic count matrix     
          }
          
          nc_dnt[d,1:max(leaves),k]=wc_dwt[d,,k] #word-count,doc d,top k on the leaves
          #recursively go up the tree to find counts on all nodes
          for(i in length(layers):1){
            for (j in 1:length(layers[[i]])){
              if (length(descendants[[layers[[i]][j]]])>0){
                nc_dnt[d,layers[[i]][j],k] = sum(nc_dnt[d,descendants[[layers[[i]][j]]],k])
              }
            }
          }
        }
        
        
      }
      
      # Initialize parameter data structures
      
      #Initialize kappa
      #kappa_pdk[p,d,k] is kappa for node p in sample d and subcommunity k
      kappa_pdk = array(0,dim=c(p,D,K))
      for(k in 1:K){
        for(a in 1:p){
          for(d in 1:D){
            kappa_pdk[a,d,k] = nc_dnt[d,descendants[[internal_nodes[a]]][1],k] - nc_dnt[d,internal_nodes[a],k]/2
          }
        }
      }
      
      #initialize phi according to a Dir(1) distribution
      #phi_dk[d,k] is the proportion of subcommunity k in sample d
      phi_dk = matrix(0,nrow=D,ncol=K)
      for(d in 1:D){
        for(k in 1:K){
          phi_dk[d,] = stats::rgamma(K,1,1)
        }
        phi_dk[d,] = phi_dk[d,]/sum(phi_dk[d,])
      }
      
      #initialize psi according to a N(0,I) distribution
      #psi_pdk[p,d,k] is the log-odds at node p in sample d and subcommunity k
      psi_pdk = array(0,dim=c(p,D,K))
      for(k in 1:K){
        for(d in 1:D){
          psi_pdk[,d,k] = chol(diag(p)) %*% matrix(stats::rnorm(p,0,1),nrow=p)
        }
      }
      
      #convert log-odds psi_pdk into probabilities theta_kda
      #theta_kda[k,d,a] is the probability in subcommunity k, sample d, and node a
      theta_kda = array(0,dim=c(K,D,A))
      for (k in 1:K){
        for (d in 1:D){
          for (a in 1:p){
            theta_kda[k,d,internal_nodes[a]] = exp(psi_pdk[a,d,k])/(1+exp(psi_pdk[a,d,k]))
          }
        }
      }
      
      #find the multinoimal distribution on the leaves implied by the theta_kda
      #beta_kdv[k,d,1:V] is the multinomial distirubtion on the V leaves in
      #subcommunity k and sample d
      beta_kdv = array(0,dim=c(K,D,V))
      for (d in 1:D){
        for (k in 1:K){
          for (leaf in leaves){ #for each leaf
            beta_kdv[k,d,leaf] = prod(theta_kda[k,d,leaf_success[[leaf]]])*prod((1-theta_kda[k,d,leaf_failures[[leaf]]]))
          }
        }
      }
      
      #initialize v_pdk, the polya-gamma auxiliary variables
      # v_pdk[p,d,k] is the PG variable associated with node p in sample d and subcommunity k
      v_pdk = array(0,dim=c(p,D,K))
      for(k in 1:K){
        for(d in 1:D){
          for(a in 1:p){
            if(nc_dnt[d,internal_nodes[a],k]<1){
              v_pdk[a,d,k] = 0
            } else {
              v_pdk[a,d,k] = pgdraw::pgdraw(nc_dnt[d,internal_nodes[a],k],psi_pdk[a,d,k])
            }
          }
        }
      }
      
      #initialize covariance matrices Sigma_ppk from the prior parameers
      #have to initialize upper and lower matrices first
      #Sigma_U_ppk[,,k] is the (diagonal) upper covariance matrix associated with subcommunity k
      Sigma_U_ppk = array(0,dim=c(p_U,p_U,K))
      for(k in 1:K){
        for(a in 1:p_U){
          Sigma_U_ppk[a,a,k] =  1/rgamma(1,shape = a_U,rate = b_U)
        }
      }
      #Initialize the graph matrix at random from the specifiedp rior
      #G_L_ppk[,,k] contains the 
      G_L_ppk = array(0,dim=c(p_L,p_L,K))
      #make data structures to store coordinates of existing/nonexisting edges
      #find the number of possible edges
      num_upptri_L = p_L*(p_L-1)/2
      #re-index the edges according to 1:num_upptri_L and find/store the coordinates
      upper_coord_L = matrix(0,nrow = num_upptri_L,ncol=2)
      upper_coord_L_C = matrix(0,nrow = num_upptri_L,ncol=2)
      #exist_ind_L_tk[t,k] records whether edge t in subcommunity k is a 0 or a 1
      exist_ind_L_tk = matrix(0,nrow = num_upptri_L,ncol = K)
      #rates_L_tk[t,k] records the birth/death rate for edge t in subcommunity k
      rates_L_tk = matrix(0,nrow = num_upptri_L,ncol = K)
      for(k in 1:K){
        coord_ct = 1
        for(i in 1:(p_L-1)){
          for(j in (i+1):p_L){
            #initialize each edge at random
            G_L_ppk[i,j,k] = sample(0:1,size = 1,prob = c(1-g_prior,g_prior))
            
            #Find the coordinates for edge coord_ct
            upper_coord_L[coord_ct,1] = i
            upper_coord_L[coord_ct,2] = j
            
            #record which edges exist
            if(G_L_ppk[i,j,k]>0){
              exist_ind_L_tk[coord_ct,k] = 1
            }
            
            #iterate counter
            coord_ct = coord_ct + 1
          }
        }
        #the diagonal is all connected
        diag(G_L_ppk[,,k]) = 1
      }
      upper_coord_L_C = upper_coord_L - 1
      #waiting_times_L_k[k] records how long the graph in subcommunity k needs to wait before an edge is born/dies
      waiting_times_L_k = rep(0,K)
      #W_L_ppk[,,k] is the lower precision matrix associated with subcommunity k
      W_L_ppk =  array(0,dim=c(p_L,p_L,K))
      for(k in 1:K){
        W_L_ppk[,,k] = f_gwish(G_L_ppk[,,k],a_L*(p_L+2),b_L*Phi_L)
      }
      #Sigma_L_ppk[,,k] is the lower covariance matrix associated with subcommunity k
      Sigma_L_ppk = array(0,dim=c(p_L,p_L,K))
      for(k in 1:K){
        Sigma_L_ppk[,,k] = solve(W_L_ppk[,,k])
      }
      #combine Sigma_k^L and Sigma_k_^U into Sigma_k
      #Sigma_ppk[,,k] is the covariance matrix associated with subcommunity k
      Sigma_ppk = array(0,dim = c(p,p,K))
      for(k in 1:K){
        Sigma_ppk[U_nodes,U_nodes,k] = Sigma_U_ppk[,,k]
        Sigma_ppk[L_nodes,L_nodes,k] = Sigma_L_ppk[,,k]
      }
      #find the inverse matrices W_ppk
      #W_ppk[,,k] is the inverse of Sigma_ppk[,,k]
      W_ppk = array(0,dim=c(p,p,K))
      for(k in 1:K){
        W_ppk[,,k] = solve(Sigma_ppk[,,k])
      }
      
      #initialize mean vectors mu_pk from a N(0,I) distribution
      #mu_pk[p,k] is the mean log-odds for node p in subcommunity k
      mu_pk = matrix(0,nrow=p,ncol=K)
      for(k in 1:K){
        mu_pk[,k] = diag(p) %*% matrix(stats::rnorm(p,0,1),nrow=p)
      }
      
      #compute Lamba Inverse
      Lambda_inv = solve(Lambda)
      
      #Pre-allocate chains
      chain_phi_dki = array(0,dim=c(D,K,iterations))
      chain_existind_L_itk = array(0,dim = c(iterations,num_upptri_L,K))
      psi_chain_k_ipd = NULL
      mu_chain_k_ip = NULL
      Sigma_chain_k_ipp = NULL
      for(k in 1:K){
        psi_chain_k_ipd[[k]] = array(0,dim=c(iterations,p,D))
        mu_chain_k_ip[[k]] = matrix(0,nrow=iterations,ncol=p) 
        Sigma_chain_k_ipp[[k]] = array(0,dim = c(iterations,p,p))
      }
      
      results = NULL
      results[[6]] = nc_dnt
      results[[5]] = chain_existind_L_itk
      results[[4]] = chain_phi_dki
      results[[3]] = psi_chain_k_ipd
      results[[2]] = mu_chain_k_ip
      results[[1]] = Sigma_chain_k_ipp
      
      cat("The burn-in period has begun. \n")
      results = LTN_Gibbs_cov_gwish_C(results, f_pg, f_gwish, Sigma_ppk, W_ppk, G_L_ppk, g_prior, upper_coord_L_C, exist_ind_L_tk, rates_L_tk, waiting_times_L_k, mu_pk, v_pdk, psi_pdk, kappa_pdk, theta_kda, beta_kdv, Lambda_inv, U_nodes_C, a_U, b_U, Phi_U, L_nodes_C, a_L, b_L, Phi_L, chain_phi_dki, psi_chain_k_ipd, mu_chain_k_ip, Sigma_chain_k_ipp, chain_existind_L_itk, nc_dnt, dt, descendants_mat_C, ta_C, docs_C, ancestors_C, internal_nodes_C, leaf_success_C, leaf_failures_C, K, p, p_U, p_L, num_upptri_L, D, V, alpha, iterations, burnin, thin)
      cat("The Gibbs Sampler has completed. \n")
      
      #unpack the results
      nc_dnt = results[[6]]
      chain_existind_L_itk = results[[5]]
      chain_phi_dki = results[[4]]
      psi_chain_k_ipd = results[[3]]
      mu_chain_k_ip = results[[2]]
      Sigma_chain_k_ipp = results[[1]]
      
      #find posterior mean for graph structure
      post_existind_L_tk = matrix(0,nrow = num_upptri_L, ncol = K)
      for(k in 1:K){
        post_existind_L_tk[,k] = apply(chain_existind_L_itk[,,k],2,mean)
      }
      #convert from the vector format to a matrix format
      post_G_L_ppk = array(0,dim = c(p_L,p_L,K))
      for(k in 1:K){
        diag(post_G_L_ppk[,,k]) = 1
        for(i in 1:num_upptri_L){
          post_G_L_ppk[upper_coord_L[i,1],upper_coord_L[i,2],k] = post_existind_L_tk[i,k]
        }
      }

      #find mean posterior psi_pdk        
      post_psi_pdk = array(0,dim=c(p,D,K))
      for(k in 1:K){
        post_psi_pd = matrix(0,nrow=p,ncol=D)
        for(d in 1:D){
          post_psi_pd[,d] = apply(psi_chain_k_ipd[[k]][1:iterations,,d],2,mean)
          # post_psi_pd[,d] = mean(psi_chain_k_ipd[[k]][1:iterations,,d])
        }
        post_psi_pdk[,,k] = post_psi_pd
      }
      
      #find theta_kdA values
      post_theta_kda = array(0,dim=c(K,D,A))
      for (k in 1:K){
        for (d in 1:D){
          for (a in 1:p){
            post_theta_kda[k,d,internal_nodes[a]] = exp(post_psi_pdk[a,d,k])/(1+exp(post_psi_pdk[a,d,k]))
          }
        }
      }
      
      #find the implicit beta-distributions
      post_beta_kdv = array(0,dim=c(K,D,V))
      for (d in 1:D){
        for (k in 1:K){
          for (leaf in leaves){ #for each leaf
            post_beta_kdv[k,d,leaf] = prod(post_theta_kda[k,d,leaf_success[[leaf]]])*prod((1-post_theta_kda[k,d,leaf_failures[[leaf]]]))
          }
        }
      }
      
      
      post_mu_pk = matrix(0,nrow=p,ncol=K)
      for(k in 1:K){
        for(a in 1:p){
          post_mu_pk[a,k] = mean(mu_chain_k_ip[[k]][,a])
        }
      }
      
      #draw theta_kdA values
      post_theta_ka = matrix(0,nrow=K,ncol=A)
      for (k in 1:K){
        for (a in 1:p){
          post_theta_ka[k,internal_nodes[a]] = exp(post_mu_pk[a,k])/(1+exp(post_mu_pk[a,k]))
        }
      }
      
      #find the implicit beta-distributions
      post_beta_kv = matrix(0,nrow=K,ncol=V)
      for (k in 1:K){
        for (leaf in leaves){ #for each leaf
          post_beta_kv[k,leaf] = prod(post_theta_ka[k,leaf_success[[leaf]]])*prod((1-post_theta_ka[k,leaf_failures[[leaf]]]))
        }
      }
      
      #psoterior Sigma means
      post_Sigma_ppk = array(0,dim=c(p,p,K))
      for(k in 1:K){
        temp = matrix(0,nrow=p,ncol=p)
        for(i in 1:iterations){
          temp = temp + Sigma_chain_k_ipp[[k]][i,,]
        }
        post_Sigma_ppk[,,k] = temp/iterations
      }
      
      #posterior phi means
      post_phi_dk = matrix(0,nrow=D,ncol=K)
      for(d in 1:D){
        for(k in 1:K){
          post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
        }
      }
      
      out = list(Mean_Post_Phi_d = post_phi_dk,
                 Mean_Post_Beta_kd = post_beta_kdv,
                 Mean_Post_Beta_k = post_beta_kv,
                 Mean_Post_G_L = post_G_L_ppk,
                 Chain_Phi = chain_phi_dki,
                 Chain_Psi = psi_chain_k_ipd,
                 Chain_Mu = mu_chain_k_ip,
                 Chain_Sigma = Sigma_chain_k_ipp,
                 phyloseq = ps
      )
      
      return(out)
      
    }
    
  }
  
  
  
  
}