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
  // priors
  // intial states
  X0 ~ normal(0,5);
  // process SD
  sigma ~ cauchy(0, 5);
  // B matrix
  diagonal(Bmat) ~ normal(0,10);
  for (i in 1:(n_spp-1)) {
    for (j in (i+1):n_spp) {
      Bmat[i, j] ~ normal(0,0.001);
      Bmat[j, i] ~ normal(0,0.001);
    }
  }
  // likelihood
  for(t in 2:n_year) {
    if(n_q == 1) {
      col(y,t) ~ normal(Bmat * col(y,t-1), sigma[1]);
    }
    else {
      col(y,t) ~ multi_normal(Bmat * col(y,t-1), diag_matrix(sigma));
    }
  }
}
