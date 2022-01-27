#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

/////////////////
// Make Multinomial Function

int sim_mult(arma::vec &prob){
  //Turn pdf vector into cdf vector
  int l = prob.size();
  arma::vec new_prob(l);
  double acc = 0;
  for (int acc_ct =0; acc_ct<l;acc_ct++){
    acc += prob[acc_ct];
    new_prob[acc_ct] = acc;
  }
  
  //Simulate multinomial from cdf vector
  NumericVector rng = runif(1);
  int ret = 0;
  for (int sim_mult_ct=0; sim_mult_ct<l; sim_mult_ct++){
    if (rng[0]<new_prob[sim_mult_ct]){
      ret = sim_mult_ct;
      break;
    }
  }
  
  return(ret);
}


// [[Rcpp::export]]
List LTN_Gibbs_C(List &results, Function f_pg,
         arma::cube &Sigma_ppk, arma::cube &W_ppk, arma::mat &mu_pk, 
         arma::cube &v_pdk, arma::cube &psi_pdk, arma::cube &kappa_pdk,
         arma::cube &theta_kda, arma::cube &beta_kdv, 
         arma::mat &Lambda_inv,
         arma::vec &gam_shape_p, arma::vec &gam_rate_p,
         arma::cube &chain_phi_dki, List &psi_chain_k_ipd, List &mu_chain_k_ip, List &Sigma_chain_k_ipp,
         arma::cube &nc_dnt, arma::mat &dt, 
         arma::mat &descendants_mat,
         List &ta, List &docs, List &ancestors,
         arma::vec &internal_nodes, List &leaf_success, List &leaf_failures,
         int &K, int &p, int &D, int &V, double & alpha,int &iterations, int &warmup, int &thin){
  
  //Initialize data structures
  
  arma::vec p_z(K);
  double max = -pow(10,20);
  arma::vec probs(K);
  arma::mat phi_dk(D,K);
  arma::mat Sigma(p,p);
  arma::mat W(p,p);
  arma::vec mu(p);
  arma::mat v(p,D);
  arma::mat psi(p,D);
  arma::mat kappa(p,D);
  NumericVector pg_draw(1);
  arma::mat v_diag(p,p);
  arma::mat psi_cov(p,p);
  arma::vec psi_mean(p);
  arma::mat psi_draw(p,1);
  arma::vec psi_bar(p);
  arma::mat mu_cov(p,p);
  arma::vec mu_mean(p);
  arma::mat draw(p,1);
  
  int it = 0;
  
  
  
  int n_it = thin*iterations + warmup;
  
  
  for(int iterate =0; iterate < n_it; iterate++){

    
    ////////////////////////////
    // Topic Assignments Loop //
    ////////////////////////////
    
    //grab phi
    //Record results
    
    for (int doc =0; doc<D; doc++){
      // int doc =0;
      NumericVector doc_temp = docs[doc];
      int N = doc_temp.size();
      
      for (int word =0; word<N;   word++){
        // int word=0;
        
        //Bookeeping section!
        //Find original topic assignment
        NumericVector old_topics = ta[doc];
        int t_old = old_topics[word];
        //Find the vocab ID of token word
        int wid = doc_temp[word];
        // //Find the ancestors of word
        arma::vec anc = ancestors[wid];
        // //Find all nodes that must be modified
        arma::vec nodes(anc.size()+1);
        for(unsigned int mk_nodes_ct =1; mk_nodes_ct<nodes.size();mk_nodes_ct++){
          nodes[mk_nodes_ct] = anc[mk_nodes_ct-1];
        }
        nodes[0] = wid;
        
        //De-increment section!
        //de-increment dt
        dt(doc,t_old) -= 1;
        //De-increment node count at old topic
        for (unsigned int nodes_ct=0; nodes_ct<nodes.size();nodes_ct++){
          nc_dnt(doc,nodes[nodes_ct],t_old) -= 1;
        }
        
        ////
        //Find probability vector for updating
        // arma::vec p_z(K);
        for(arma::uword p_z_ct =0; p_z_ct<p_z.size();p_z_ct++){
          p_z[p_z_ct] = log((dt(doc,p_z_ct) + alpha)) + log(beta_kdv(p_z_ct,doc,wid));
        }
        
        max = -pow(10,20);
        for(unsigned int max_ct=0; max_ct<p_z.size(); max_ct++){
          if (p_z[max_ct]>max){
            max = p_z[max_ct];
          }
        }
        
        probs = exp(p_z - max)/sum(exp(p_z-max));
        // arma::vec probs(K);
        // for(int k=0; k<K; k++){
        //   probs[k] = 0.5;
        // }
        
        //Draw multinomial
        int t_new = sim_mult(probs); //Need to readjust types of vector
        
        // Re-increment section!
        
        //Update topic assignments
        old_topics[word] = t_new;
        ta[doc] = old_topics;
        
        // Re-increment dt
        dt(doc,t_new) += 1;
        // Re-increment node count at new topic
        for (arma::uword nodes_ct=0; nodes_ct<nodes.size();nodes_ct++){
          nc_dnt(doc,nodes[nodes_ct],t_new) += 1;
        }
        
        
      }
      
      
    }
    
    //////////////////
    // Update Kappa //
    //////////////////
    for(int k=0; k<K; k++){
      // int k =0;
      for(int a=0; a<p; a++){
        // int a = 0;
        for(int d=0; d<D; d++){
          // int d = 0 ;
          int parent = internal_nodes[a];
          int child = descendants_mat(parent,0);
          kappa_pdk(a,d,k) = nc_dnt(d,child,k) - nc_dnt(d,parent,k)/2;
        }
      }
    }
    
    
    
    //////////////
    // LTN Loop //
    /////////////
    for(int k=0; k<K; k++){
      // int k = 0;
      Sigma = Sigma_ppk.slice(k);
      W = W_ppk.slice(k);
      mu = mu_pk.col(k);
      v = v_pdk.slice(k);
      psi = psi_pdk.slice(k);
      kappa = kappa_pdk.slice(k);
      
      
      //////////////
      // Update v //
      //////////////
      
      for(int d=0; d<D; d++){
        // int d = 0;
        for(int a=0; a<p; a++){
          // int a = 0;
          int b = nc_dnt(d,internal_nodes[a],k);
          // int b = 0;
          if (b<1){
            v(a,d) = 0 ;
          } else if (b < 30){
            double c = psi(a,d);
            pg_draw = f_pg(b,c);
            v(a,d) = pg_draw[0];
          } else {
            double c = psi(a,d);
            double pg_mean = b/(2*c)*tanh(c/2);
            double pg_sd = sqrt(b/(4*pow(c,3))*(1/pow(cosh(c/2),2))*(sinh(c)-c));
            
            pg_draw = Rcpp::rnorm(1,pg_mean,pg_sd);
            v(a,d) = pg_draw[0];
          }
        }
      }
      
      
      ////////////////
      // Update Psi //
      ////////////////
      
      arma::vec avg = W * mu;
      
      for(int d=0; d<D; d++){
        // int d = 0;
        v_diag = v_diag.zeros();
        for(int v_diag_ct=0; v_diag_ct<p; v_diag_ct++){
          v_diag(v_diag_ct,v_diag_ct) = v(v_diag_ct,d);
        }
        
        psi_cov = W + v_diag;
        psi_cov = arma::symmatu(psi_cov);
        psi_cov = arma::inv(psi_cov); //Do in two steps for memory?
        psi_mean = psi_cov * (avg + kappa.col(d));
        psi_draw = Rcpp::rnorm(p,0,1);
        
        //Do in two steps for memory?
        psi_cov = arma::symmatu(psi_cov);
        psi_cov = chol(psi_cov);
        psi.col(d) = psi_cov*psi_draw + psi_mean;
      }
      psi_pdk.slice(k) = psi;
      
      //Find the average value of psi for each node a
      for(int a=0; a<p; a++){
        double sum = 0;
        for(int d=0; d < D; d++){
          sum += psi(a,d);
        }
        psi_bar[a] = sum/D;
      }
      
      ///////////////
      // Update Mu //
      ///////////////
      mu_cov = Lambda_inv + D*W;
      mu_cov = arma::symmatu(mu_cov);
      mu_cov = arma::inv(mu_cov); //do in two steps for memory
      mu_mean = ((mu_cov) * W ) * (psi_bar) * D;
      
      //do some rearranging for memory reasons
      mu_cov = arma::symmatu(mu_cov);
      mu_cov = chol(mu_cov);
      
      draw = Rcpp::rnorm(p,0,1);
      mu = mu_cov*draw + mu_mean;
      mu_pk.col(k) = mu;
      
      ///////////////////////
      // Update Sigma and W//
      ///////////////////////
      Sigma = Sigma.zeros();
      for(int sig_a=0; sig_a<p; sig_a++){
        
        double Sigma_sum = 0;
        for(int sig_d=0; sig_d<D; sig_d++){
          Sigma_sum += pow(psi(sig_a,sig_d) - mu[sig_a],2);
        }
        
        arma::vec tau_draw = Rcpp::rgamma(1,gam_shape_p[sig_a] + D/2, 2/(2*gam_rate_p[sig_a] + Sigma_sum));
        Sigma(sig_a,sig_a) = 1/tau_draw[0];
      }
      Sigma = arma::symmatu(Sigma);
      Sigma_ppk.slice(k) = Sigma;
      W = arma::inv(Sigma);
      W = arma::symmatu(W);
      W_ppk.slice(k) = W;
      
    }
    
    /////////////////////////
    // Convert psi to beta //
    /////////////////////////
    
    //Convert psi to theta
    for(int k=0; k<K; k++){
      for(int d=0; d<D; d++){
        for(int a=0; a<p; a++){
          theta_kda(k,d,internal_nodes[a]) = exp(psi_pdk(a,d,k))/(1+exp(psi_pdk(a,d,k)));
        }
      }
    }
    
    //Convert theta to beta
    for(int k=0; k<K; k++){
      // int k = 0;
      for(int d=0; d<D; d++){
        // int d = 0;
        for(int v=0; v<V; v++){
          // int v = 1;
          arma::vec success_ind = leaf_success[v];
          double num_suc = success_ind.size();
          double prod = 1;
          for(int suc_ct=0; suc_ct<num_suc; suc_ct++){
            if (success_ind[suc_ct] < 0){
              
            } else {
              prod = prod*theta_kda(k,d,success_ind[suc_ct]);
            }
          }
          
          arma::vec fail_ind = leaf_failures[v];
          double num_fail = fail_ind.size();
          for(int fail_ct=0; fail_ct<num_fail; fail_ct++){
            if (fail_ind[fail_ct] < 0){
              
            } else {
              prod = prod * (1-theta_kda(k,d,fail_ind[fail_ct]));
            }
          }
          
          beta_kdv(k,d,v) = prod;
          
        }
      }
    }
    
    if (iterate < warmup){
      if(iterate  % 100 ==0){
        if(iterate > 0){
          Rcout << "The burn-in iterate is " << iterate << " of " << warmup << ".\n";
        }
      }
    }
    
    if(iterate >= warmup){
      
      if(iterate <= warmup){
        Rcout << "The burn-in period has completed. \n";
        Rcout << "The sampling has begun. \n";
      }
      
      if(iterate % thin == 0 ){
        /////////////////////
        // Record Results! //
        /////////////////////
        
        /////////
        // Phi //
        /////////
        // arma::mat phi_dk = chain_phi_dki.slice(it);
        phi_dk = chain_phi_dki.slice(it);
        for(int d=0; d<D; d++){
          //Find the sum for the denominator
          double sum = 0;
          for(int k = 0; k<K; k++){
            sum += dt(d,k);
          }
          
          //Find the esitmated value of phi
          for(int k = 0; k<K; k++){
            phi_dk(d,k) = (dt(d,k) + alpha)/(sum + K*alpha);
          }
        }
        chain_phi_dki.slice(it) = phi_dk;
        
        //////////////
        // LTN Loop //
        /////////////
        for(int k=0; k<K; k++){
          // int k = 0;
          Sigma = Sigma_ppk.slice(k);
          W = W_ppk.slice(k);
          mu = mu_pk.col(k);
          v = v_pdk.slice(k);
          psi = psi_pdk.slice(k);
          kappa = kappa_pdk.slice(k);
          
          
          ////////////////
          // Record Psi //
          ////////////////
          arma::cube temp_psi = psi_chain_k_ipd[k];
          temp_psi.row(it) = psi;
          psi_chain_k_ipd[k] = temp_psi;
          
          ///////////////
          // Record Mu //
          ///////////////
          arma::mat temp_mu = mu_chain_k_ip[k];
          temp_mu.row(it) = trans(mu);
          mu_chain_k_ip[k] = temp_mu;
          
          //////////////////
          // Record Sigma //
          //////////////////
          arma::cube temp_sig = Sigma_chain_k_ipp[k];
          temp_sig.row(it) = Sigma;
          Sigma_chain_k_ipp[k] = temp_sig;
          
        }
        
        if(it % 10 ==0){
          if(it>0){
            Rcout << "The Gibbs Sampler iterate is " << it << " of " << iterations << ".\n";
          }
        }
        
        it += 1;

        
      }
    }
    
    
    
    
    
  }
  
  results[4] = nc_dnt;
  results[3] = chain_phi_dki;
  results[2] = psi_chain_k_ipd;
  results[1] = mu_chain_k_ip;
  results[0] = Sigma_chain_k_ipp;

  
  return(results);
}

// [[Rcpp::export]]
List LTN_Gibbs_cov_C(List &results, Function f_pg, Function f_iwish,
                 arma::cube &Sigma_ppk, arma::cube &W_ppk, arma::mat &mu_pk, 
                 arma::cube &v_pdk, arma::cube &psi_pdk, arma::cube &kappa_pdk,
                 arma::cube &theta_kda, arma::cube &beta_kdv, 
                 arma::mat &Lambda_inv,
                 arma::vec &U_nodes, double &a_U, double &b_U, arma::mat &Phi_U,
                 arma::vec &L_nodes, double &a_L, double &b_L, arma::mat &Phi_L,
                 arma::cube &chain_phi_dki, List &psi_chain_k_ipd, List &mu_chain_k_ip, List &Sigma_chain_k_ipp,
                 arma::cube &nc_dnt, arma::mat &dt, 
                 arma::mat &descendants_mat,
                 List &ta, List &docs, List &ancestors,
                 arma::vec &internal_nodes, List &leaf_success, List &leaf_failures,
                 int &K, int &p, int &p_U, int &p_L, int &D, int &V, double & alpha,int &iterations, int &warmup, int &thin){
  
  //Initialize data structures
  
  arma::vec p_z(K);
  double max = -pow(10,20);
  arma::vec probs(K);
  arma::mat phi_dk(D,K);
  arma::mat Sigma(p,p);
  arma::mat W(p,p);
  arma::vec mu(p);
  arma::mat v(p,D);
  arma::mat psi(p,D);
  arma::mat kappa(p,D);
  NumericVector pg_draw(1);
  arma::mat v_diag(p,p);
  arma::mat psi_cov(p,p);
  arma::vec psi_mean(p);
  arma::mat psi_draw(p,1);
  arma::vec psi_bar(p);
  arma::mat mu_cov(p,p);
  arma::vec mu_mean(p);
  arma::mat draw(p,1);
  
  //New data structures
  NumericMatrix Sigma_U(p_U,p_U);
  arma::mat Sigma_U_mat(p_U,p_U);
  arma::mat psi_U_pd(p_U,D);
  arma::vec mu_U(p_U);
  NumericMatrix Sigma_L(p_L,p_L);
  arma::mat Sigma_L_mat(p_L,p_L);
  arma::mat psi_L_pd(p_L,D);
  arma::vec mu_L(p_L);

  
  int it = 0;
  
  
  
  int n_it = thin*iterations + warmup;
  
  
  for(int iterate =0; iterate < n_it; iterate++){
    
    
    ////////////////////////////
    // Topic Assignments Loop //
    ////////////////////////////
    
    //grab phi
    //Record results
    
    for (int doc =0; doc<D; doc++){
      // int doc =0;
      NumericVector doc_temp = docs[doc];
      int N = doc_temp.size();
      
      for (int word =0; word<N;   word++){
        // int word=0;
        
        //Bookeeping section!
        //Find original topic assignment
        NumericVector old_topics = ta[doc];
        int t_old = old_topics[word];
        //Find the vocab ID of token word
        int wid = doc_temp[word];
        // //Find the ancestors of word
        arma::vec anc = ancestors[wid];
        // //Find all nodes that must be modified
        arma::vec nodes(anc.size()+1);
        for(unsigned int mk_nodes_ct =1; mk_nodes_ct<nodes.size();mk_nodes_ct++){
          nodes[mk_nodes_ct] = anc[mk_nodes_ct-1];
        }
        nodes[0] = wid;
        
        //De-increment section!
        //de-increment dt
        dt(doc,t_old) -= 1;
        //De-increment node count at old topic
        for (unsigned int nodes_ct=0; nodes_ct<nodes.size();nodes_ct++){
          nc_dnt(doc,nodes[nodes_ct],t_old) -= 1;
        }
        
        ////
        //Find probability vector for updating
        // arma::vec p_z(K);
        for(arma::uword p_z_ct =0; p_z_ct<p_z.size();p_z_ct++){
          p_z[p_z_ct] = log((dt(doc,p_z_ct) + alpha)) + log(beta_kdv(p_z_ct,doc,wid));
        }
        
        max = -pow(10,20);
        for(unsigned int max_ct=0; max_ct<p_z.size(); max_ct++){
          if (p_z[max_ct]>max){
            max = p_z[max_ct];
          }
        }
        
        probs = exp(p_z - max)/sum(exp(p_z-max));
        // arma::vec probs(K);
        // for(int k=0; k<K; k++){
        //   probs[k] = 0.5;
        // }
        
        //Draw multinomial
        int t_new = sim_mult(probs); //Need to readjust types of vector
        
        // Re-increment section!
        
        //Update topic assignments
        old_topics[word] = t_new;
        ta[doc] = old_topics;
        
        // Re-increment dt
        dt(doc,t_new) += 1;
        // Re-increment node count at new topic
        for (arma::uword nodes_ct=0; nodes_ct<nodes.size();nodes_ct++){
          nc_dnt(doc,nodes[nodes_ct],t_new) += 1;
        }
        
        
      }
      
      
    }
    
    //////////////////
    // Update Kappa //
    //////////////////
    for(int k=0; k<K; k++){
      // int k =0;
      for(int a=0; a<p; a++){
        // int a = 0;
        for(int d=0; d<D; d++){
          // int d = 0 ;
          int parent = internal_nodes[a];
          int child = descendants_mat(parent,0);
          kappa_pdk(a,d,k) = nc_dnt(d,child,k) - nc_dnt(d,parent,k)/2;
        }
      }
    }
    
    
    
    //////////////
    // LTN Loop //
    /////////////
    for(int k=0; k<K; k++){
      // int k = 0;
      Sigma = Sigma_ppk.slice(k);
      W = W_ppk.slice(k);
      mu = mu_pk.col(k);
      v = v_pdk.slice(k);
      psi = psi_pdk.slice(k);
      kappa = kappa_pdk.slice(k);
      
      
      //////////////
      // Update v //
      //////////////
      
      for(int d=0; d<D; d++){
        // int d = 0;
        for(int a=0; a<p; a++){
          // int a = 0;
          int b = nc_dnt(d,internal_nodes[a],k);
          // int b = 0;
          if (b<1){
            v(a,d) = 0 ;
          } else if (b < 30){
            double c = psi(a,d);
            pg_draw = f_pg(b,c);
            v(a,d) = pg_draw[0];
          } else {
            double c = psi(a,d);
            double pg_mean = b/(2*c)*tanh(c/2);
            double pg_sd = sqrt(b/(4*pow(c,3))*(1/pow(cosh(c/2),2))*(sinh(c)-c));
            
            pg_draw = Rcpp::rnorm(1,pg_mean,pg_sd);
            v(a,d) = pg_draw[0];
          }
        }
      }
      
      
      ////////////////
      // Update Psi //
      ////////////////
      
      arma::vec avg = W * mu;
      
      for(int d=0; d<D; d++){
        // int d = 0;
        v_diag = v_diag.zeros();
        for(int v_diag_ct=0; v_diag_ct<p; v_diag_ct++){
          v_diag(v_diag_ct,v_diag_ct) = v(v_diag_ct,d);
        }
        
        psi_cov = W + v_diag;
        psi_cov = arma::symmatu(psi_cov);
        psi_cov = arma::inv(psi_cov); //Do in two steps for memory?
        psi_mean = psi_cov * (avg + kappa.col(d));
        psi_draw = Rcpp::rnorm(p,0,1);
        
        //Do in two steps for memory?
        psi_cov = arma::symmatu(psi_cov);
        psi_cov = chol(psi_cov);
        psi.col(d) = psi_cov*psi_draw + psi_mean;
      }
      psi_pdk.slice(k) = psi;
      
      //Find the average value of psi for each node a
      for(int a=0; a<p; a++){
        double sum = 0;
        for(int d=0; d < D; d++){
          sum += psi(a,d);
        }
        psi_bar[a] = sum/D;
      }
      
      ///////////////
      // Update Mu //
      ///////////////
      mu_cov = Lambda_inv + D*W;
      mu_cov = arma::symmatu(mu_cov);
      mu_cov = arma::inv(mu_cov); //do in two steps for memory
      mu_mean = ((mu_cov) * W ) * (psi_bar) * D;
      
      //do some rearranging for memory reasons
      mu_cov = arma::symmatu(mu_cov);
      mu_cov = chol(mu_cov);
      
      draw = Rcpp::rnorm(p,0,1);
      mu = mu_cov*draw + mu_mean;
      mu_pk.col(k) = mu;
      
      ///////////////////////
      // Update Sigma and W//
      ///////////////////////
      
      /////////////
      // Sigma_U //
      /////////////
      
      //Find psi_U and mu_U
      for(int sig_U_node = 0; sig_U_node<p_U; sig_U_node++){
        psi_U_pd.row(sig_U_node) = psi.row(U_nodes[sig_U_node]);
        mu_U[sig_U_node] = mu[U_nodes[sig_U_node]];
      }
      
      //Find matrix for posterior
      Sigma_U_mat = Sigma_U_mat.zeros();
      for(int sig_d = 0; sig_d<D; sig_d++){
        Sigma_U_mat += (psi_U_pd.col(sig_d) - mu_U) * trans(psi_U_pd.col(sig_d) - mu_U);
      }
      
      //Draw matrix 
      //NB not an arma --- have to switch if we wish to parallelize this, but its not a bottleneck
      Sigma_U = f_iwish(D + a_U*(p_U+2),Sigma_U_mat + b_U*Phi_U);
      
      /////////////
      // Sigma_L //
      /////////////
      
      //Find psi_L and mu_L
      for(int sig_L_node = 0; sig_L_node<p_L; sig_L_node++){
        psi_L_pd.row(sig_L_node) = psi.row(L_nodes[sig_L_node]);
        mu_L[sig_L_node] = mu[L_nodes[sig_L_node]];
      }
      
      //Find matrix for posterior
      Sigma_L_mat = Sigma_L_mat.zeros();
      for(int sig_d = 0; sig_d<D; sig_d++){
        Sigma_L_mat += (psi_L_pd.col(sig_d) - mu_L) * trans(psi_L_pd.col(sig_d) - mu_L);
      }
      
      //Draw matrix
      //NB not an arma --- have to switch if we wish to parallelize this, but its not a bottleneck
      Sigma_L = f_iwish(D + a_L*(p_L+2),Sigma_L_mat + b_L*Phi_L);
      
      ///////////
      // Sigma //
      ///////////
      
      Sigma = Sigma.zeros();
      
      //put Sigma_U into Sigma
      for(int sig_U_i = 0; sig_U_i<p_U; sig_U_i++){
        for(int sig_U_ii = 0; sig_U_ii<p_U; sig_U_ii++){
          Sigma(U_nodes[sig_U_i],U_nodes[sig_U_ii]) = Sigma_U(sig_U_i,sig_U_ii);
        }
      }
      //put Sigma_L into Sigma
      for(int sig_L_i = 0; sig_L_i<p_L; sig_L_i++){
        for(int sig_L_ii = 0; sig_L_ii<p_L; sig_L_ii++){
          Sigma(L_nodes[sig_L_i],L_nodes[sig_L_ii]) = Sigma_L(sig_L_i,sig_L_ii);
        }
      }
      
      
      Sigma = arma::symmatu(Sigma);
      Sigma_ppk.slice(k) = Sigma;
      W = arma::inv(Sigma);
      W = arma::symmatu(W);
      W_ppk.slice(k) = W;
      
    }
    
    /////////////////////////
    // Convert psi to beta //
    /////////////////////////
    
    //Convert psi to theta
    for(int k=0; k<K; k++){
      for(int d=0; d<D; d++){
        for(int a=0; a<p; a++){
          theta_kda(k,d,internal_nodes[a]) = exp(psi_pdk(a,d,k))/(1+exp(psi_pdk(a,d,k)));
        }
      }
    }
    
    //Convert theta to beta
    for(int k=0; k<K; k++){
      // int k = 0;
      for(int d=0; d<D; d++){
        // int d = 0;
        for(int v=0; v<V; v++){
          // int v = 1;
          arma::vec success_ind = leaf_success[v];
          double num_suc = success_ind.size();
          double prod = 1;
          for(int suc_ct=0; suc_ct<num_suc; suc_ct++){
            if (success_ind[suc_ct] < 0){
              
            } else {
              prod = prod*theta_kda(k,d,success_ind[suc_ct]);
            }
          }
          
          arma::vec fail_ind = leaf_failures[v];
          double num_fail = fail_ind.size();
          for(int fail_ct=0; fail_ct<num_fail; fail_ct++){
            if (fail_ind[fail_ct] < 0){
              
            } else {
              prod = prod * (1-theta_kda(k,d,fail_ind[fail_ct]));
            }
          }
          
          beta_kdv(k,d,v) = prod;
          
        }
      }
    }
    
    if (iterate < warmup){
      if(iterate  % 100 ==0){
        if(iterate > 0){
          Rcout << "The burn-in iterate is " << iterate << " of " << warmup << ".\n";
        }
      }
    }
    
    if(iterate >= warmup){
      
      if(iterate <= warmup){
        Rcout << "The burn-in period has completed. \n";
        Rcout << "The sampling has begun. \n";
      }
      
      if(iterate % thin == 0 ){
        /////////////////////
        // Record Results! //
        /////////////////////
        
        /////////
        // Phi //
        /////////
        // arma::mat phi_dk = chain_phi_dki.slice(it);
        phi_dk = chain_phi_dki.slice(it);
        for(int d=0; d<D; d++){
          //Find the sum for the denominator
          double sum = 0;
          for(int k = 0; k<K; k++){
            sum += dt(d,k);
          }
          
          //Find the esitmated value of phi
          for(int k = 0; k<K; k++){
            phi_dk(d,k) = (dt(d,k) + alpha)/(sum + K*alpha);
          }
        }
        chain_phi_dki.slice(it) = phi_dk;
        
        //////////////
        // LTN Loop //
        /////////////
        for(int k=0; k<K; k++){
          // int k = 0;
          Sigma = Sigma_ppk.slice(k);
          W = W_ppk.slice(k);
          mu = mu_pk.col(k);
          v = v_pdk.slice(k);
          psi = psi_pdk.slice(k);
          kappa = kappa_pdk.slice(k);
          
          
          ////////////////
          // Record Psi //
          ////////////////
          arma::cube temp_psi = psi_chain_k_ipd[k];
          temp_psi.row(it) = psi;
          psi_chain_k_ipd[k] = temp_psi;
          
          ///////////////
          // Record Mu //
          ///////////////
          arma::mat temp_mu = mu_chain_k_ip[k];
          temp_mu.row(it) = trans(mu);
          mu_chain_k_ip[k] = temp_mu;
          
          //////////////////
          // Record Sigma //
          //////////////////
          arma::cube temp_sig = Sigma_chain_k_ipp[k];
          temp_sig.row(it) = Sigma;
          Sigma_chain_k_ipp[k] = temp_sig;
          
        }
        
        if(it % 10 ==0){
          if(it>0){
            Rcout << "The Gibbs Sampler iterate is " << it << " of " << iterations << ".\n";
          }
        }
        
        it += 1;
        
        
      }
    }
    
    
    
    
    
  }
  
  results[4] = nc_dnt;
  results[3] = chain_phi_dki;
  results[2] = psi_chain_k_ipd;
  results[1] = mu_chain_k_ip;
  results[0] = Sigma_chain_k_ipp;
  
  
  return(results);
}

// [[Rcpp::export]]
List LTN_Gibbs_Perp_C(List &results, Function f_pg,
         arma::cube &Sigma_ppk, arma::cube &W_ppk, arma::mat &mu_pk, 
         arma::cube &v_pdk, arma::cube &psi_pdk, arma::cube &kappa_pdk,
         arma::cube &theta_kda, arma::cube &beta_kdv, 
         arma::cube &chain_phi_dki, List &psi_chain_k_ipd,
         arma::cube &nc_dnt, arma::mat &dt, 
         arma::mat &descendants_mat,
         List &ta, List &docs, List &ancestors,
         arma::vec &internal_nodes, List &leaf_success, List &leaf_failures,
         int &K, int &p, int &D, int &V, double & alpha,int &iterations, int &warmup, int &thin){
  
  //Initialize data structures
  
  arma::vec p_z(K);
  double max = -pow(10,20);
  arma::vec probs(K);
  arma::mat phi_dk(D,K);
  arma::mat Sigma(p,p);
  arma::mat W(p,p);
  arma::vec mu(p);
  arma::mat v(p,D);
  arma::mat psi(p,D);
  arma::mat kappa(p,D);
  NumericVector pg_draw(1);
  arma::mat v_diag(p,p);
  arma::mat psi_cov(p,p);
  arma::vec psi_mean(p);
  arma::mat psi_draw(p,1);
  arma::vec psi_bar(p);
  arma::mat mu_cov(p,p);
  arma::vec mu_mean(p);
  arma::mat draw(p,1);
  
  int it = 0;
  
  
  
  int n_it = thin*iterations + warmup;
  
  
  for(int iterate =0; iterate < n_it; iterate++){

    
    ////////////////////////////
    // Topic Assignments Loop //
    ////////////////////////////
    
    //grab phi
    //Record results
    
    for (int doc =0; doc<D; doc++){
      // int doc =0;
      NumericVector doc_temp = docs[doc];
      int N = doc_temp.size();
      
      for (int word =0; word<N;   word++){
        // int word=0;
        
        //Bookeeping section!
        //Find original topic assignment
        NumericVector old_topics = ta[doc];
        int t_old = old_topics[word];
        //Find the vocab ID of token word
        int wid = doc_temp[word];
        // //Find the ancestors of word
        arma::vec anc = ancestors[wid];
        // //Find all nodes that must be modified
        arma::vec nodes(anc.size()+1);
        for(unsigned int mk_nodes_ct =1; mk_nodes_ct<nodes.size();mk_nodes_ct++){
          nodes[mk_nodes_ct] = anc[mk_nodes_ct-1];
        }
        nodes[0] = wid;
        
        //De-increment section!
        //de-increment dt
        dt(doc,t_old) -= 1;
        //De-increment node count at old topic
        for (unsigned int nodes_ct=0; nodes_ct<nodes.size();nodes_ct++){
          nc_dnt(doc,nodes[nodes_ct],t_old) -= 1;
        }
        
        ////
        //Find probability vector for updating
        // arma::vec p_z(K);
        for(arma::uword p_z_ct =0; p_z_ct<p_z.size();p_z_ct++){
          p_z[p_z_ct] = log((dt(doc,p_z_ct) + alpha)) + log(beta_kdv(p_z_ct,doc,wid));
        }
        
        max = -pow(10,20);
        for(unsigned int max_ct=0; max_ct<p_z.size(); max_ct++){
          if (p_z[max_ct]>max){
            max = p_z[max_ct];
          }
        }
        
        probs = exp(p_z - max)/sum(exp(p_z-max));
        // arma::vec probs(K);
        // for(int k=0; k<K; k++){
        //   probs[k] = 0.5;
        // }
        
        //Draw multinomial
        int t_new = sim_mult(probs); //Need to readjust types of vector
        
        // Re-increment section!
        
        //Update topic assignments
        old_topics[word] = t_new;
        ta[doc] = old_topics;
        
        // Re-increment dt
        dt(doc,t_new) += 1;
        // Re-increment node count at new topic
        for (arma::uword nodes_ct=0; nodes_ct<nodes.size();nodes_ct++){
          nc_dnt(doc,nodes[nodes_ct],t_new) += 1;
        }
        
        
      }
      
      
    }
    
    //////////////////
    // Update Kappa //
    //////////////////
    for(int k=0; k<K; k++){
      // int k =0;
      for(int a=0; a<p; a++){
        // int a = 0;
        for(int d=0; d<D; d++){
          // int d = 0 ;
          int parent = internal_nodes[a];
          int child = descendants_mat(parent,0);
          kappa_pdk(a,d,k) = nc_dnt(d,child,k) - nc_dnt(d,parent,k)/2;
        }
      }
    }
    
    
    
    //////////////
    // LTN Loop //
    /////////////
    for(int k=0; k<K; k++){
      // int k = 0;
      Sigma = Sigma_ppk.slice(k);
      W = W_ppk.slice(k);
      mu = mu_pk.col(k);
      v = v_pdk.slice(k);
      psi = psi_pdk.slice(k);
      kappa = kappa_pdk.slice(k);
      
      
      //////////////
      // Update v //
      //////////////
      
      for(int d=0; d<D; d++){
        // int d = 0;
        for(int a=0; a<p; a++){
          // int a = 0;
          int b = nc_dnt(d,internal_nodes[a],k);
          // int b = 0;
          if (b<1){
            v(a,d) = 0 ;
          } else if (b < 30){
            double c = psi(a,d);
            pg_draw = f_pg(b,c);
            v(a,d) = pg_draw[0];
          } else {
            double c = psi(a,d);
            double pg_mean = b/(2*c)*tanh(c/2);
            double pg_sd = sqrt(b/(4*pow(c,3))*(1/pow(cosh(c/2),2))*(sinh(c)-c));
            
            pg_draw = Rcpp::rnorm(1,pg_mean,pg_sd);
            v(a,d) = pg_draw[0];
          }
        }
      }
      
      
      ////////////////
      // Update Psi //
      ////////////////
      
      arma::vec avg = W * mu;
      
      for(int d=0; d<D; d++){
        // int d = 0;
        v_diag = v_diag.zeros();
        for(int v_diag_ct=0; v_diag_ct<p; v_diag_ct++){
          v_diag(v_diag_ct,v_diag_ct) = v(v_diag_ct,d);
        }
        
        psi_cov = W + v_diag;
        psi_cov = arma::symmatu(psi_cov);
        psi_cov = arma::inv(psi_cov); //Do in two steps for memory?
        psi_mean = psi_cov * (avg + kappa.col(d));
        psi_draw = Rcpp::rnorm(p,0,1);
        
        //Do in two steps for memory?
        psi_cov = arma::symmatu(psi_cov);
        psi_cov = chol(psi_cov);
        psi.col(d) = psi_cov*psi_draw + psi_mean;
      }
      psi_pdk.slice(k) = psi;
      
      //Find the average value of psi for each node a
      for(int a=0; a<p; a++){
        double sum = 0;
        for(int d=0; d < D; d++){
          sum += psi(a,d);
        }
        psi_bar[a] = sum/D;
      }
      
    }
    
    /////////////////////////
    // Convert psi to beta //
    /////////////////////////
    
    //Convert psi to theta
    for(int k=0; k<K; k++){
      for(int d=0; d<D; d++){
        for(int a=0; a<p; a++){
          theta_kda(k,d,internal_nodes[a]) = exp(psi_pdk(a,d,k))/(1+exp(psi_pdk(a,d,k)));
        }
      }
    }
    
    //Convert theta to beta
    for(int k=0; k<K; k++){
      // int k = 0;
      for(int d=0; d<D; d++){
        // int d = 0;
        for(int v=0; v<V; v++){
          // int v = 1;
          arma::vec success_ind = leaf_success[v];
          double num_suc = success_ind.size();
          double prod = 1;
          for(int suc_ct=0; suc_ct<num_suc; suc_ct++){
            if (success_ind[suc_ct] < 0){
              
            } else {
              prod = prod*theta_kda(k,d,success_ind[suc_ct]);
            }
          }
          
          arma::vec fail_ind = leaf_failures[v];
          double num_fail = fail_ind.size();
          for(int fail_ct=0; fail_ct<num_fail; fail_ct++){
            if (fail_ind[fail_ct] < 0){
              
            } else {
              prod = prod * (1-theta_kda(k,d,fail_ind[fail_ct]));
            }
          }
          
          beta_kdv(k,d,v) = prod;
          
        }
      }
    }
    
    if (iterate < warmup){
      if(iterate % 100 ==0){
        if(iterate > 0){
          Rcout << "The burn-in iterate is " << iterate << " of " << warmup << ".\n";
        }
      }
    }
    
    if(iterate >= warmup){
      
      if(iterate <= warmup){
        Rcout << "The burn-in period has completed. \n";
        Rcout << "The sampling has begun. \n";
      }
      
      if(iterate % thin == 0 ){
        /////////////////////
        // Record Results! //
        /////////////////////
        
        /////////
        // Phi //
        /////////
        // arma::mat phi_dk = chain_phi_dki.slice(it);
        phi_dk = chain_phi_dki.slice(it);
        for(int d=0; d<D; d++){
          //Find the sum for the denominator
          double sum = 0;
          for(int k = 0; k<K; k++){
            sum += dt(d,k);
          }
          
          //Find the esitmated value of phi
          for(int k = 0; k<K; k++){
            phi_dk(d,k) = (dt(d,k) + alpha)/(sum + K*alpha);
          }
        }
        chain_phi_dki.slice(it) = phi_dk;
        
        //////////////
        // LTN Loop //
        /////////////
        for(int k=0; k<K; k++){
          // int k = 0;
          psi = psi_pdk.slice(k);

          
          ////////////////
          // Record Psi //
          ////////////////
          arma::cube temp_psi = psi_chain_k_ipd[k];
          temp_psi.row(it) = psi;
          psi_chain_k_ipd[k] = temp_psi;
          
        }
        
        if(it % 10 ==0){
          if(it >0){
            Rcout << "The Gibbs Sampler iterate is " << it << " of " << iterations << ".\n";
          }
        }
        
        it += 1;
        
      }
    }
    
    
    
    
    
  }
  
  results[2] = nc_dnt;
  results[1] = chain_phi_dki;
  results[0] = psi_chain_k_ipd;

  
  return(results);
}