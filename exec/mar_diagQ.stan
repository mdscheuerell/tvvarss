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
  matrix[n_spp,n_year] y;
}
parameters {
  matrix<lower=-1,upper=1>[n_spp,n_spp] Bmat; // B matrix
  vector<lower=0>[n_spp] sigma;               // proc SD
  vector[n_spp] X0;                           // initial states
}
transformed parameters {

}
model {
  X0 ~ normal(0,5);
  sigma ~ cauchy(0, 5); // half-Cauchy due to constraint
  for(i 2:n_year) {
    if(n_q > 1) {
      y[i] ~ normal()
    }
  }
}
