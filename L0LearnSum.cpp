// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 

using namespace Rcpp;
using namespace arma;

// global variables to easily pass through functions
double M; // adaptive tuning
int L0_norm; double L1_norm, L2_norm; // norms
int totit, numit; // iterations

///////////////////////
/// Data Processing ///
///////////////////////
// return beta_out on original scale
// [[Rcpp::export]]
void rescale(mat & beta_out, vec & scaling) {
  int p = beta_out.n_rows;
  int n_par = beta_out.n_cols;
  
  // rescale in place
  for (int j = 0; j < p; j++) {
    for (int k = 0; k < n_par; k++) {
      beta_out(j,k) = beta_out(j,k) * scaling(j);
    }
  }
}

//////////////////////////
/// Coordinate Descent ///
//////////////////////////
// [[Rcpp::export]]
double L0LearnSum(vec & beta0, const mat & R, const vec & r, const vec & ix,
                  double lambda0, double lambda1, double lambda2, 
                  const bool adap=true, int maxit=1e2, double btol=1e-4) {
  int p = beta0.size(); 
  int j = 0;
  
  double threshold = sqrt((2*lambda0) / (1+2*lambda2));
  
  double dbeta = 0.0, max_dbeta = 0.0;
  double beta_j = 0.0, beta_tilde = 0.0, z = 0.0;
  
  vec R_beta = R * beta0; // R beta
  int L0 = 0; double L1 = 0.0, L2 = 0.0;
  
  double conv = 0; numit = maxit;
  
  for (int t = 0; t < maxit; t++) {
    max_dbeta = 0.0; 
    L0 = 0; L1 = 0.0; L2 = 0.0;
    
    for (int j_ix = 0; j_ix < p; j_ix++) {
      j = ix(j_ix);
      
      beta_j = beta0(j); beta0(j) = 0.0;
      
      // coordinate update
      beta_tilde = r(j) - (R_beta(j) - beta_j); // R=X'X, r=X'Y
      z = (std::abs(beta_tilde) - lambda1) / (1+2*lambda2);
      if (z > threshold) beta0(j) = std::copysign(z, beta_tilde);
      
      // hard-code fix for exploding beta
      if (std::abs(beta0(j)) > 10.0) beta0(j) = r(j);

      dbeta = beta0(j) - beta_j;
      
      // update if necessary
      if (std::abs(dbeta) > 1e-20) {
        // max change 
        max_dbeta = std::max(max_dbeta, std::abs(dbeta));
        
        // products with beta
        R_beta += R.col(j) * dbeta; 
      }
      
      // norms
      if (std::abs(beta0(j)) > 1e-20) L0 += 1;
      L1 += std::abs(beta0(j));
      L2 += std::pow(beta0(j), 2.0);
    }
    
    // check stability of active sets (if relevant)
    if (max_dbeta <= btol) {
      numit = t+1; 
      conv = 1;
      break;
    }
  }
  
  // update for norms
  L0_norm += L0;
  L1_norm += L1;
  L2_norm += L2;
  
  // adaptive tuning step after convergence (skip otherwise)
  if (adap & (conv == 1)) {
    double m0 = 0.0;
    for (int j = 0; j < p; j++) {
      if (std::abs(beta0(j)) <= 1e-20) {
        beta_tilde = r(j) - R_beta(j);
        m0 = std::pow(std::abs(beta_tilde) - lambda1, 2.0) / (2 * (1 + 2*lambda2));
        if (m0 > M) M = m0;
      } 
    }
  }
  
  return conv;
}

// [[Rcpp::export]]
double L0LearnSum_active(vec & beta0, const mat & R, const vec & r, const vec & ix,
                         double lambda0, double lambda1, double lambda2, 
                         int max_active=1, const bool adap=true, int maxit=1e2) {
  int p = beta0.size(); 
  int j = 0;
  
  bool active_changed = false;
  int active_set_number = max_active;
  
  double threshold = sqrt((2*lambda0) / (1+2*lambda2));
  
  double dbeta = 0.0;
  // double max_dbeta = 0.0;
  double beta_j = 0.0, beta_tilde = 0.0, z = 0.0;
  
  vec R_beta = R * beta0; // R beta
  int L0 = 0; double L1 = 0.0, L2 = 0.0;
  
  double conv = 0; numit = maxit;
  
  for (int t = 0; t < maxit; t++) {
    // max_dbeta = 0.0; 
    L0 = 0; L1 = 0.0; L2 = 0.0;
    active_changed = false;
    
    for (int j_ix = 0; j_ix < p; j_ix++) {
      j = ix(j_ix);
      
      // skip if active set hasn't converged yet and beta_j = 0
      if ((active_set_number < max_active) & (std::abs(beta0(j)) < 1e-20)) continue;
      
      beta_j = beta0(j); beta0(j) = 0.0;
      
      // coordinate update
      beta_tilde = r(j) - (R_beta(j) - beta_j); // R=X'X, r=X'Y
      z = (std::abs(beta_tilde) - lambda1) / (1+2*lambda2);
      if (z > threshold) beta0(j) = std::copysign(z, beta_tilde);
      
      // hard-code fix for exploding beta
      if (std::abs(beta0(j)) > 10.0) beta0(j) = r(j);
      
      dbeta = beta0(j) - beta_j;
      
      // update if necessary
      if (std::abs(dbeta) > 1e-20) {
        // max change 
        // max_dbeta = std::max(max_dbeta, std::abs(dbeta));
        
        // products with beta
        R_beta += R.col(j) * dbeta; 
      } 
      
      // active set has changed (either newly zeroed or selected)
      if (((std::abs(beta_j) > 1e-20) & (std::abs(beta0(j)) < 1e-20)) |
          ((std::abs(beta_j) < 1e-20) & (std::abs(beta0(j)) > 1e-20))) {
        active_changed = true;
      }
      
      // norms
      if (std::abs(beta0(j)) > 1e-20) L0 += 1;
      L1 += std::abs(beta0(j));
      L2 += std::pow(beta0(j), 2.0);
    }
    
    // check if active set if stabilized
    if (!active_changed) {
      active_set_number++;
    } else {
      active_set_number = 0;
    }
    
    // check stability of active sets + convergence
    if (active_set_number == max_active) {
      numit = t+1; 
      conv = 1;
      break;
    }
  }
  
  // update for norms
  L0_norm += L0;
  L1_norm += L1;
  L2_norm += L2;
  
  // adaptive tuning step after convergence (skip otherwise)
  if ((adap) & (conv == 1)) {
    double m0 = 0.0;
    for (int j = 0; j < p; j++) {
      if (std::abs(beta0(j)) <= 1e-20) {
        beta_tilde = r(j) - R_beta(j);
        m0 = std::pow(std::abs(beta_tilde) - lambda1, 2.0) / (2 * (1 + 2*lambda2));
        if (m0 > M) M = m0;
      } 
    }
  }
  // Rcout << numit << '\n';
  
  return conv;
}

//////////////
/// Blocks ///
//////////////
// [[Rcpp::export]]
double L0LearnSum_block(vec & beta_vec, const List & R_list, const vec & r_full, 
                        const List & ix_list, const vec & start, const vec & stop,
                        double lambda0, double lambda1, double lambda2, 
                        const bool adap=true, int maxit=1e2, double btol=1e-4) {
  int n_blocks = R_list.size();
  int p = beta_vec.size();
  double conv_full = 0.0, conv0 = 0.0;
  totit = 0;
  
  // loop over blocks
  for (int k=0; k<n_blocks; k++) {
    vec beta_block = beta_vec.subvec(start(k), stop(k));
    conv0 = L0LearnSum(beta_block, R_list(k), r_full.subvec(start(k), stop(k)), ix_list(k),
                       lambda0, lambda1, lambda2, 
                       adap, maxit, btol);
    
    beta_vec.subvec(start(k), stop(k)) = beta_block;
    
    conv_full += conv0;

    totit = std::max(totit, numit);
    
  }
  
  return conv_full / n_blocks;
}

// [[Rcpp::export]]
double L0LearnSum_active_block(vec & beta_vec, const List & R_list, const vec & r_full, 
                               const List & ix_list, const vec & start, const vec & stop,
                               double lambda0, double lambda1, double lambda2, 
                               int max_active=1, const bool adap=true, int maxit=1e2, double btol=1e-4) {
  int n_blocks = R_list.size();
  int p = beta_vec.size();
  double conv_full = 0.0, conv0 = 0.0;
  totit = 0;
  
  // loop over blocks
  for (int k=0; k<n_blocks; k++) {
    vec beta_block = beta_vec.subvec(start(k), stop(k));
    conv0 = L0LearnSum_active(beta_block, R_list(k), r_full.subvec(start(k), stop(k)), ix_list(k),
                              lambda0, lambda1, lambda2, 
                              max_active, adap, maxit);
    
    beta_vec.subvec(start(k), stop(k)) = beta_block;
    
    conv_full += conv0;
    
    totit = std::max(totit, numit);
    
  }
  
  return conv_full / n_blocks;
}


/////////////////////
/// Adaptive Grid ///
/////////////////////
// par_grid = [lambda0_start, lambda1, lambda2, alpha]
// [[Rcpp::export]]
mat L0LearnSum_auto(mat & beta_out, const List & R_list, const vec & r_full, 
                    const List & ix_list, vec & start, vec & stop,
                    mat par_grid, int n_lambda0, vec scaling=0,
                    int maxit=1e2, double btol=1e-4) {
  int p = r_full.size();
  int n_par = par_grid.n_rows * n_lambda0;
  
  int n_blocks = R_list.size();
  
  mat par(n_par, 9); // [lambda0, lambda1, lambda2, M, L0, L1, L2, conv, totit]
  vec beta_full(p);
  
  double converged = 0.0;
  
  double lambda0=0.0; double lambda1=0.0; double lambda2=0.0; double alpha=0.0; 
  int grid_ix = 0;
  int progress = 1; int ten_pct = n_par / 10;
  
  M = 0.0;
  for (int par_ix = 0; par_ix < par_grid.n_rows; par_ix++) {
    lambda0 = par_grid(par_ix, 0);
    lambda1 = par_grid(par_ix, 1);
    lambda2 = par_grid(par_ix, 2);
    alpha = par_grid(par_ix, 3);
    
    // restart beta at every new lambda1/lambda2
    beta_full.fill(0);
    
    for (int lambda0_ix = 0; lambda0_ix < n_lambda0; lambda0_ix++) {
      if ((grid_ix+1) % ten_pct == 0) {
        Rcout << progress << "/10..";
        progress++;
      }
      
      M = 0.0; L0_norm = 0; L1_norm = 0.0; L2_norm = 0.0;
      converged = L0LearnSum_block(beta_full, R_list, r_full, ix_list, start, stop, 
                                   lambda0, lambda1, lambda2,  
                                   true, maxit, btol);
      
      par(grid_ix, 0) = lambda0; par(grid_ix, 1) = lambda1; par(grid_ix, 2) = lambda2;
      par(grid_ix, 3) = M; par(grid_ix, 4) = L0_norm;
      par(grid_ix, 5) = L1_norm; par(grid_ix, 6) = L2_norm;
      par(grid_ix, 7) = converged; par(grid_ix, 8) = totit;
      
      beta_out.col(grid_ix) = beta_full;
      
      if (M >= lambda0 | std::abs(M) < 1e-20) {
        lambda0 = alpha * lambda0;
      } else {
        lambda0 = alpha * M;
      }
      
      grid_ix++;
    }
  }
  
  // rescale beta
  if (scaling.size() == p) {
    Rcout << "rescaling";
    rescale(beta_out, scaling);
  }
  
  Rcout << "\n";
  
  return par;
}

// [[Rcpp::export]]
mat L0LearnSum_active_auto(mat & beta_out, const List & R_list, const vec & r_full, 
                           const List & ix_list, vec & start, vec & stop,
                           mat par_grid, int n_lambda0, vec scaling=0,
                           int max_active=1, int maxit=1e2, double btol=1e-4) {
  int p = r_full.size();
  int n_par = par_grid.n_rows * n_lambda0;
  
  int n_blocks = R_list.size();
  
  mat par(n_par, 9); // [lambda0, lambda1, lambda2, M, L0, L1, L2, conv, totit]
  vec beta_full(p);
  
  double converged = 0.0;
  
  double lambda0=0.0; double lambda1=0.0; double lambda2=0.0; double alpha=0.0; 
  int grid_ix = 0;
  int progress = 1; int ten_pct = n_par / 10;
  
  M = 0.0;
  for (int par_ix = 0; par_ix < par_grid.n_rows; par_ix++) {
    lambda0 = par_grid(par_ix, 0);
    lambda1 = par_grid(par_ix, 1);
    lambda2 = par_grid(par_ix, 2);
    alpha = par_grid(par_ix, 3);
    
    // restart beta at every new lambda1/lambda2
    beta_full.fill(0);
    
    for (int lambda0_ix = 0; lambda0_ix < n_lambda0; lambda0_ix++) {
      if ((grid_ix+1) % ten_pct == 0) {
        Rcout << progress << "/10..";
        progress++;
      }
      
      M = 0.0; L0_norm = 0; L1_norm = 0.0; L2_norm = 0.0;
      converged = L0LearnSum_active_block(beta_full, R_list, r_full, ix_list, start, stop, 
                                          lambda0, lambda1, lambda2,  
                                          max_active, true, maxit);
      
      par(grid_ix, 0) = lambda0; par(grid_ix, 1) = lambda1; par(grid_ix, 2) = lambda2;
      par(grid_ix, 3) = M; par(grid_ix, 4) = L0_norm;
      par(grid_ix, 5) = L1_norm; par(grid_ix, 6) = L2_norm;
      par(grid_ix, 7) = converged; par(grid_ix, 8) = totit;
      
      beta_out.col(grid_ix) = beta_full;
      
      if (M >= lambda0 | std::abs(M) < 1e-20) {
        lambda0 = alpha * lambda0;
      } else {
        lambda0 = alpha * M;
      }
      
      grid_ix++;
    }
  }
  
  // rescale beta
  if (scaling.size() == p) {
    Rcout << "rescaling";
    rescale(beta_out, scaling);
  }
  
  Rcout << "\n";
  
  return par;
}


////////////////
// Fixed Grid //
////////////////
// par_grid = [lambda0, lambda1, lambda2, alpha]
// lambda0 decreasing sequencing
// [[Rcpp::export]]
mat L0LearnSum_grid(mat & beta_out, const List & R_list, const vec & r_full, 
                    const List & ix_list, vec & start, vec & stop,
                    mat par_grid, vec scaling=0,
                    int maxit=1e2, double btol=1e-4) {
  int p = r_full.size();
  int n_par = par_grid.n_rows;
  
  int n_blocks = R_list.size();
  
  mat par(n_par, 9); // [lambda0, lambda1, lambda2, M, L0, L1, L2, conv, totit]
  vec beta_full(p);
  
  double converged = 0.0;
  
  double s=0.0; double lambda0=0.0; double lambda1=0.0; double lambda2=0.0;
  int progress = 1; int ten_pct = n_par / 10;
  
  M = 0.0;
  for (int par_ix = 0; par_ix < n_par; par_ix++) {
    if ((par_ix+1) % ten_pct == 0) {
      Rcout << progress << "/10..";
      progress++;
    }
    
    // restart beta with a new lambda1/lambda2
    if (par_grid(par_ix, 0) > lambda0) beta_full.fill(0);
    
    lambda0 = par_grid(par_ix, 0);
    lambda1 = par_grid(par_ix, 1);
    lambda2 = par_grid(par_ix, 2);
    
    L0_norm = 0; L1_norm = 0.0; L2_norm = 0.0;
    converged = L0LearnSum_block(beta_full, R_list, r_full, ix_list, start, stop, 
                                 lambda0, lambda1, lambda2, 
                                 false, maxit, btol);
    
    par(par_ix, 0) = lambda0; par(par_ix, 1) = lambda1; par(par_ix, 2) = lambda2;
    par(par_ix, 3) = M; par(par_ix, 4) = L0_norm;
    par(par_ix, 5) = L1_norm; par(par_ix, 6) = L2_norm;
    par(par_ix, 7) = converged; par(par_ix, 8) = totit;
    
    beta_out.col(par_ix) = beta_full;
  }
  
  // rescale beta
  if (scaling.size() == p) {
    Rcout << "rescaling";
    rescale(beta_out, scaling);
  } 
  
  Rcout << "\n";
  
  return par;
}

mat L0LearnSum_active_grid(mat & beta_out, const List & R_list, const vec & r_full, 
                           const List & ix_list, vec & start, vec & stop,
                           mat par_grid, vec scaling=0,
                           int max_active=1, int maxit=1e2, double btol=1e-4) {
  int p = r_full.size();
  int n_par = par_grid.n_rows;
  
  int n_blocks = R_list.size();
  
  mat par(n_par, 9); // [lambda0, lambda1, lambda2, M, L0, L1, L2, conv, totit]
  vec beta_full(p);
  
  double converged = 0.0;
  
  double s=0.0; double lambda0=0.0; double lambda1=0.0; double lambda2=0.0;
  int progress = 1; int ten_pct = n_par / 10;
  
  M = 0.0;
  for (int par_ix = 0; par_ix < n_par; par_ix++) {
    if ((par_ix+1) % ten_pct == 0) {
      Rcout << progress << "/10..";
      progress++;
    }
    
    // restart beta with a new lambda1/lambda2
    if (par_grid(par_ix, 0) > lambda0) beta_full.fill(0);
    
    lambda0 = par_grid(par_ix, 0);
    lambda1 = par_grid(par_ix, 1);
    lambda2 = par_grid(par_ix, 2);
    
    L0_norm = 0; L1_norm = 0.0; L2_norm = 0.0;
    converged = L0LearnSum_active_block(beta_full, R_list, r_full, ix_list, start, stop, 
                                        lambda0, lambda1, lambda2, 
                                        max_active, false, maxit, btol);
    
    par(par_ix, 0) = lambda0; par(par_ix, 1) = lambda1; par(par_ix, 2) = lambda2;
    par(par_ix, 3) = M; par(par_ix, 4) = L0_norm;
    par(par_ix, 5) = L1_norm; par(par_ix, 6) = L2_norm;
    par(par_ix, 7) = converged; par(par_ix, 8) = totit;
    
    beta_out.col(par_ix) = beta_full;
  }
  
  // rescale beta
  if (scaling.size() == p) {
    Rcout << "rescaling";
    rescale(beta_out, scaling);
  } 
  
  Rcout << "\n";
  
  return par;
}


