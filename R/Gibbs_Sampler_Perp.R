#' Runs a modified LTN-LDA Gibbs sampler to estimate perplexity on a test set.  
#'
#' This function takes the output of an LTNLDA function ran on the training set, an edge matrix specifying a tree (which should be the same tree as for the training set), and a matrix specifying counts of sequencing reads per sample on the test set. It runs a modified collapsed blocked Gibbs sampler for LTNLDA for these inputs to estimate the phi_{d} and psi_{p,d,k} on the test set, and returns a list containing a perplexity estimate on the test set as well as  posterior mean estimates and Markov Chains for the parameters.  For a detailed example see the vignette "Perplexity".
#'
#' @param model The output of an LTNLDA model trained on a training set.  
#' @param ps A phyloseq object corresponding to the test set containing an otu_table() and a phy_tree() with an edge matrix.  That is, otu_table(ps) and phy_tree(ps)$edge both exist.
#' @param iterations The number of iterations to record.  Default value is 1000.
#' @param burnin The number of burnin iterations to run before recording values.  The default value is 10000.
#' @param thin The amount by which we thin the chain.  A value of X means that 1 every X values is recorded after the burnin.  The default value is 10.
#' @param alpha A double specifying the prior on the subcommunity-sample proportions.  The default value is 1.
#' @return A list with 5 entries. Perplexity is a double which is the estimated perplexity on the test set.  Mean_Post_Phi_d contains the posterior mean estimate for the subcommunity-sample distributions phi_d on each sample in the test set.  Mean_Post_Psi_pdk contains the posterior mean estimates for the log-odds psi_{p,d,k} for the test set. Chain_Phi contains the Markov Chains for the phi_d. Chain_Psi contains the Markov Chains for the psi_{p,d,k}.
#' @examples
#' perp = LTNLDA_Perplexity(model, test_ps)
#' @references
#' ADD PAPER INFORMATION ONCE WE KNOW IT
#' @export

LTNLDA_Perplexity = function(model, ps,
                iterations = 1000, burnin = 10000, thin = 10, alpha = 1){
  
  #check inputs to see if they have been specified or not
  if(is.null(model)){
    cat("Model is not specified. \n")
  } else{
    if(is.null(ps)){
      cat("Phyloseq object is not specified. \n")
    } else{
      
      cat("Processing data. \n")
            
            #extract structures from phyloseq
            tree.edge = phyloseq::phy_tree(ps)$edge
            test_dtm = phyloseq::otu_table(ps)
            #load K from model 
            K = dim(model$Final_Iterate_Counts)[3]
            
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
            

            
            
            #######################################################
            # Convert dtm count data to list of vectors of tokens #
            #######################################################

            #find number of samples
            D = ncol(test_dtm)

            #Initialize an extra count matrix
            extra_dtm = matrix(0,nrow=V,ncol=D)

            #This code partitions the count matrix about in half - into a test matrix and an extra matrix -
            #For document completion
            #We fit the model on the test_dtm and find perplexity on extra_dtm

            #if there are less than 3 sequencing reads, it stays in the test set
            #This doesn't matter because the count is so low
            #If it is more than three, we split about in half
            #This is a janky way to handle even and odd counts
            for(d in 1:D){
              for(v in 1:V){
                total = test_dtm[v,d]
                if( total > 3){
                  extra_dtm[v,d] = sample( (round(total/2) - 1):(round(total/2) + 1),1)
                  test_dtm[v,d] = test_dtm[v,d]-extra_dtm[v,d]
                }

              }
            }


            #convert test_dtm count data to list of vectors of tokens
            docs = NULL
            docs[[D+1]] = rep(0,10000)
            docs_C = NULL
            for (d in 1:D){ #for each document
              doc = NULL

              for(v in 1:V){
                doc = c(doc,rep(v,test_dtm[v,d]))
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


            ###################################
            # Load posterior means from model #
            ###################################
            post_Sigma_ppk = model$Mean_Post_Sigma_ppk
            post_mu_pk = model$Mean_Post_Mu_k
            #Initialize global parameters to these mean values



            #initialize covariance matrices Sigma_ppk from the posterior means
            #Sigma_ppk[,,k] is the covariance matrix associated with subcommunity k
            Sigma_ppk = array(0,dim=c(p,p,K))
            for(k in 1:K){
                Sigma_ppk[,,k] = post_Sigma_ppk[,,k]
            }
            #find the inverse matrices W_ppk
            #W_ppk[,,k] is the inverse of Sigma_ppk[,,k]
            W_ppk = array(0,dim=c(p,p,K))
            for(k in 1:K){
              for(a in 1:p){
                W_ppk[a,a,k] = 1/Sigma_ppk[a,a,k]
              }
            }

            #initialize mean vectors mu_pk from the posterior means
            #mu_pk[p,k] is the mean log-odds for node p in subcommunity k
            mu_pk = matrix(0,nrow=p,ncol=K)
            for(k in 1:K){
              mu_pk[,k] = post_mu_pk[,k]
            }
            

            #Pre-allocate chains
            chain_phi_dki = array(0,dim=c(D,K,iterations))
            psi_chain_k_ipd = NULL
            for(k in 1:K){
              psi_chain_k_ipd[[k]] = array(0,dim=c(iterations,p,D))
            }

            chains = NULL
            chains[[3]] = nc_dnt
            chains[[2]] = chain_phi_dki
            chains[[1]] =   psi_chain_k_ipd

            cat("The burn-in period has begun. \n")
            chains = LTN_Gibbs_Perp_C(chains, f_pg, Sigma_ppk, W_ppk, mu_pk, v_pdk, psi_pdk, kappa_pdk, theta_kda, beta_kdv,  chain_phi_dki, psi_chain_k_ipd, nc_dnt, dt, descendants_mat_C, ta_C, docs_C, ancestors_C, internal_nodes_C, leaf_success_C, leaf_failures_C, K, p, D, V, alpha, iterations, burnin, thin)
            cat("The Gibbs Sampler has completed. \n")
            
            nc_dnt = chains[[3]]
            chain_phi_dki = chains[[2]]
            psi_chain_k_ipd = chains[[1]]

            #Find the phi_d on the test set
            test_post_phi_dk = matrix(0,nrow=D,ncol=K)
            for(d in 1:D){
              for(k in 1:K){
                test_post_phi_dk[d,k] = mean(chain_phi_dki[d,k,1:iterations])
              }
            }
            
            #find the posterior mean estimates for psi_pdk onthe test set
            test_post_psi_pdk = array(0,dim=c(p,D,K))
            for(k in 1:K){
              test_post_psi_pd = matrix(0,nrow=p,ncol=D)
              for(d in 1:D){
                test_post_psi_pd[,d] = apply(psi_chain_k_ipd[[k]][1:iterations,,d],2,mean)
              }
              test_post_psi_pdk[,,k] = test_post_psi_pd
            }

            #find the value of theta and beta for each iteration on the test set
            post_theta_kdai = array(0,dim=c(K,D,A,iterations))
            post_beta_kdvi = array(0,dim=c(K,D,V,iterations))
            for(i in 1:iterations){
              #convert psi to theta
              for (k in 1:K){
                for (d in 1:D){
                  for (a in 1:p){
                    post_theta_kdai[k,d,internal_nodes[a],i] = exp(psi_chain_k_ipd[[k]][i,a,d])/(1+exp(psi_chain_k_ipd[[k]][i,a,d]))
                  }
                }
              }
              #convert theta to beta
              for (d in 1:D){
                for (k in 1:K){
                  for (leaf in leaves){ #for each leaf
                    post_beta_kdvi[k,d,leaf,i] = prod(post_theta_kdai[k,d,leaf_success[[leaf]],i])*prod((1-post_theta_kdai[k,d,leaf_failures[[leaf]],i]))
                  }
                }
              }
            }


            #######################
            # Evaluate Perplexity #
            #######################

            #convert extra_dtm into extra_docs
            extra_docs = NULL
            extra_docs[[D+1]] = rep(0,10000)
            for (d in 1:D){ #for each document
              doc = NULL

              for(v in 1:V){
                doc = c(doc,rep(v,extra_dtm[v,d]))
              }

              doc = doc[sample(1:length(doc),replace = FALSE)]
              extra_docs[[d]] = doc
            }
            extra_docs[[D+1]] = NULL

            #Estimate perplexity on the extra set
            doc_comp_prob = rep(0,D)

            for (d in 1:D){
              total = 0
              for (i in 1:iterations){
                subtotal =0
                for (w in 1:length(extra_docs[[d]])){
                  wid = extra_docs[[d]][w]

                  subtotal = subtotal + log(sum(chain_phi_dki[d,,i]*post_beta_kdvi[,d,wid,i]))

                }
                total = subtotal+total
              }
              doc_comp_prob[d] = total/(iterations)
            }

            den = 0
            for (d in 1:length(extra_docs)){
              den = den + length(extra_docs[[d]])
            }
            perp = exp(-sum(doc_comp_prob)/den)
            
            out = list(Perplexity = perp,
                       Mean_Post_Phi_d = test_post_phi_dk,
                       Mean_Post_Psi_pdk = test_post_psi_pdk,
                       Chain_Phi = chain_phi_dki,
                       Chain_Psi = psi_chain_k_ipd)
            
            return(out)
        }
  }
  
  
  
}