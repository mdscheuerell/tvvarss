functions {
  // mapB maps b onto [lo,up]
  // mapB(b, 0, 1) = inverse_logit(b)
  real mapB(real b, real lo, real up) {
    real b_star;
    b_star = (up - lo) / (1 + exp(-b)) + lo;
    return b_star;
  }
}
data {
  // scalars
  int<lower=0> n_year;                      // # years
  int<lower=0> n_site;                      // # sites
  int<lower=0> n_spp;                       // # species
  // int<lower=0> n_process;                   // # total processes
  int<lower=0> n_q;                         // # proc variances = max(shared_q)
  int<lower=0> n_q_c;                       // # covariances in Q
  // vectors
  int<lower=0> b_diag[n_spp*n_spp];         // is this on diag of B: Y/N = 2/1
  int<lower=0> b_indices[n_spp*n_spp];      // indices for B values
  int<lower=0> row_indices[n_spp*n_spp];    // row indices for B values
  int<lower=0> col_indices[n_spp*n_spp];    // col indices for B values
  int<lower=0> spp_indices_pos[n_pos];      // species IDs
  int<lower=0> year_indices_pos[n_pos];     // year IDs
  real y[n_pos];                            // real-valued data
  // matrices
  int<lower=0> shared_q[n_spp,n_process+1]; // which sites/spp share proc variances
  int<lower=0> shared_r[n_spp,n_site+1];    // which sites/spp share obs variances
  int<lower=0> shared_u[n_spp,n_process+1]; // which sites/spp share trends
  matrix[n_process,n_spp] x0;               // initial states
  matrix[4,2] b_limits;                     // lo/up constraints on B elements
}
parameters {
  matrix<lower=-1,upper=1>[n_spp,n_spp] Bmat; // B matrix
  corr_matrix[n_spp] Omega; // correlation matrix
  vector<lower=0>[n_spp] sigma; // scales  vector
  X0[n_spp];                            // initial states
}
transformed parameters {
  if(n_q_c > 0) {
    cov_matrix[n_spp] Sigma; // covariance matrix
    // for (m in 1:n_spp) {
    //   Sigma[m, m] = sigma[m] * sigma[m] * Omega[m, m];
    // }
    diagonal(Sigma) = diag_matrix(sigma) * diag_matrix(sigma) * Omega
    for (i in 1:(n_spp-1)) {
      for (j in (i+1):n_spp) {
        Sigma[i, j] = sigma[i] * sigma[j] * Omega[i, j];
        Sigma[j, i] = Sigma[i, j];
      }
    }
  }
}
      else {
        Sigma[i, j] = 0
      }
model {
  if(n_q_c > 0) {
    Omega ~ lkj_corr(2.0); // regularize to unit correlation
    sigma ~ cauchy(0, 5); // half-Cauchy due to constraint
    for(t in 2:n_year) {

    }
  }
}
