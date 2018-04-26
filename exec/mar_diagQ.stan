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
  int<lower=0> n_spp;                       // # species
  // int<lower=0> n_process;                   // # total processes
  int<lower=0> n_q;                         // # proc variances = max(shared_q)
  // matrices
  matrix[n_spp,n_year] y;
}
parameters {
  matrix<lower=-1,upper=1>[n_spp,n_spp] Bmat; // B matrix
  vector<lower=0>[n_spp] sigma;               // proc SD
  vector[n_spp] X0;                           // initial states
}
model {
  X0 ~ normal(0,5);
  sigma ~ cauchy(0, 5); // half-Cauchy due to constraint
  for(t 2:n_year) {
    if(n_q == 1) {
      col(y,t) ~ normal(col(y,t-1),sigma[1])
    }
    else {
      col(y,t) ~ multi_normal(Bmat * col(y,t-1), diag_matrix(sigma))
    }
  }
}
