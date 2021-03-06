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
  int<lower=0> n_process;                   // # total processes
  int<lower=0> n_q;                         // # proc variances = max(shared_q)
  int<lower=0> n_r;                         // # obs variances = max(shared_r)
  int<lower=0> n_u;                         // # trends = max(shared_u)
  int<lower=0> n_pos;                       // # non-NA values in the data
  int<lower=0> fit_dynamicB;                // fit dynamic B? Y/N = 1/0
  int<lower=0> est_trend;                   // bias in RW? Y/N = 1/0
  int<lower=0> demean;                      // non-zero mean inproc? Y/N = 1/0
  int family;                               // distn form of obs data
  // vectors
  int<lower=0> process[n_site+1];           // map of sites:processes
  int<lower=0> b_diag[n_spp*n_spp];         // is this on diag of B: Y/N = 2/1
  int<lower=0> b_indices[n_spp*n_spp];      // indices for B values
  int<lower=0> row_indices[n_spp*n_spp];    // row indices for B values
  int<lower=0> col_indices[n_spp*n_spp];    // col indices for B values
  int<lower=0> spp_indices_pos[n_pos];      // species IDs
  int<lower=0> site_indices_pos[n_pos];     // site:proc IDs
  int<lower=0> year_indices_pos[n_pos];     // year IDs
  real y[n_pos];                            // real-valued data
  int y_int[n_pos];                         // integer data
  // matrices
  int<lower=0> shared_q[n_spp,n_process+1]; // which sites/spp share proc variances
  int<lower=0> shared_r[n_spp,n_site+1];    // which sites/spp share obs variances
  int<lower=0> shared_u[n_spp,n_process+1]; // which sites/spp share trends
  matrix[n_process,n_spp] x0;               // initial states
  matrix[4,2] b_limits;                     // lo/up constraints on B elements
}
transformed data {
  int<lower=0> n_spp2;
  n_spp2 = n_spp * n_spp;
}
parameters {
  // vector<lower=-1,upper=1>[n_spp2] vecBdev[n_year]; // elements accessed [n_year,n_spp]
  vector<lower=-1,upper=1>[n_spp2] vecBdev; // elements accessed [n_year,n_spp]
  // real<lower=0> sigma_rw_pars[2];                   // sds for random walk
  matrix[n_year,n_spp] x[n_process];                // unobserved states
  real<lower=0> resid_process_sd[n_q];              // SD of proc errors
  real<lower=0> obs_sd[n_r];                        // SD of obs errors
  real u[n_u];                                      // biases/trends in RW
}
transformed parameters {
  // vector<lower=0>[n_spp2] sigma_rw;
  matrix<lower=0>[n_spp, n_process] resid_process_mat;
  matrix<lower=0>[n_spp, n_site] obs_mat;
  vector<lower=-20,upper=20>[n_spp2] vecB[n_year];
  matrix[n_spp, n_process] u_mat;
  matrix[n_spp,n_spp] B[(n_year-1)]; // B matrix, accessed as n_year, n_spp, n_spp
  matrix[n_year,n_spp] pred[n_process]; // predicted unobserved states

  // for(i in 1:n_spp2) {
  //   sigma_rw[i] = sigma_rw_pars[b_diag[i]];
  // }
  for(i in 1:n_spp) {
    for(j in 1:n_process) {
      resid_process_mat[i,j] = resid_process_sd[shared_q[i,j]];
      u_mat[i,j] = u[shared_u[i,j]];
    }
    for(j in 1:n_site) {
      obs_mat[i,j] = obs_sd[shared_r[i,j]];
    }
  }

  for(s in 1:n_process) {
    pred[s,1,] = x[s,1,]; // states for first year
  }
  // for(i in 1:n_spp2) {
    // vecB[1,i] = vecBdev[1,i]; // first time step, prior in model {}
  // }
  vecB[1] = vecBdev;
  for(t in 2:n_year) {
    // fill in B matrix, shared across sites

    // estimated interactions
    for(i in 1:n_spp2) {
      // map b_ij to appropriate prior range
      B[t-1,row_indices[i],col_indices[i]] = mapB(vecB[t-1,i], b_limits[b_indices[i],1], b_limits[b_indices[i],2]);
      // random walk in b elements
      // if user wants constant b matrix, vecB[t] = vecB[t-1]
      // vecB[t,i] = vecB[t-1,i] + fit_dynamicB * vecBdev[t-1,i];
      vecB[t,i] = vecB[t-1,i];
    }

    // do projection to calculate predicted values. modify code depending on whether
    // predictions should be demeaned before projected, and whether or not trend included.
    for(s in 1:n_process) {
     if(est_trend == 0) {
      if(demean==0) {pred[s,t,] = x[s,t-1,] * B[t-1,,];}
      if(demean==1) {pred[s,t,] = (x[s,t-1,] - u_mat[,s]') * B[t-1,,];}
     }
     if(est_trend == 1) {
      if(demean==0) {pred[s,t,] = x[s,t-1,] * B[t-1,,] + u_mat[,s]';}
      if(demean==1) {pred[s,t,] = (x[s,t-1,] - u_mat[,s]') * B[t-1,,] + u_mat[,s]';}
     }
    }
  }
}
model {
  // sigma_rw_pars[1] ~ student_t(5,0,1);
  // sigma_rw_pars[2] ~ student_t(5,0,1);
  // for(t in 1:n_year) {
  //   vecBdev[t] ~ normal(0, 1); // vectorized random in B // removed (0, sigma_rw)
  // }
  for(i in 1:n_spp2) {
      vecBdev[i] ~ normal(0, 10);
  }
  // process model
  for(site in 1:n_process) {
    for(spp in 1:n_spp) {
      // time = 1
      x[site,1,spp] ~ normal(x0[site,spp],1);
      // time = 2+
      for(t in 2:n_year) {
        x[site,t,spp] ~ normal(pred[site,t,spp], resid_process_mat[spp,site]);
      }
    }
  }
  // prior on process standard deviations
  for(i in 1:n_q) {
    resid_process_sd[i] ~ student_t(5,0,1);
  }
  // prior on obs standard deviations
  for(i in 1:n_r) {
    obs_sd[i] ~ student_t(5,0,1);
  }
  // prior on biases/trends
  for(i in 1:n_u) {
    u[i] ~ normal(0,1);
  }
  // observation model
  // gaussian likelihood for now
  for(i in 1:n_pos) {
    if(family==1) y[i] ~ normal(x[site_indices_pos[i],year_indices_pos[i],spp_indices_pos[i]], obs_mat[spp_indices_pos[i],site_indices_pos[i]]);
    if(family==2) y_int[i] ~ bernoulli_logit(x[site_indices_pos[i],year_indices_pos[i],spp_indices_pos[i]]);
    if(family==3) y_int[i] ~ poisson_log(x[site_indices_pos[i],year_indices_pos[i],spp_indices_pos[i]]);
    if(family==4) y[i] ~ gamma(obs_mat[spp_indices_pos[i],site_indices_pos[i]], obs_mat[spp_indices_pos[i],site_indices_pos[i]] ./ x[site_indices_pos[i],year_indices_pos[i],spp_indices_pos[i]]);
    if(family==5) y[i] ~ lognormal(x[site_indices_pos[i],year_indices_pos[i],spp_indices_pos[i]], obs_mat[spp_indices_pos[i],site_indices_pos[i]]);
  }

}
generated quantities {
  vector[n_spp] B_fix;
  vector[n_pos] log_lik;
  // diag(B_0)
  for(i in 1:n_spp) {
    B_fix[i] = B[1,i,i];
  }
  // for use in loo() package
  if(family==1) for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | x[site_indices_pos[n],year_indices_pos[n],spp_indices_pos[n]], obs_mat[spp_indices_pos[n],site_indices_pos[n]]);
  if(family==2) for (n in 1:n_pos) log_lik[n] = bernoulli_lpmf(y_int[n] | inv_logit(x[site_indices_pos[n],year_indices_pos[n],spp_indices_pos[n]]) );
  if(family==3) for (n in 1:n_pos) log_lik[n] = poisson_lpmf(y_int[n] | exp(x[site_indices_pos[n],year_indices_pos[n],spp_indices_pos[n]]) );
  if(family==4) for (n in 1:n_pos) log_lik[n] = gamma_lpdf(y[n] | obs_mat[spp_indices_pos[n],site_indices_pos[n]], obs_mat[spp_indices_pos[n],site_indices_pos[n]] ./ x[site_indices_pos[n],year_indices_pos[n],spp_indices_pos[n]]);
  if(family==5) for (n in 1:n_pos) log_lik[n] = lognormal_lpdf(y[n] | x[site_indices_pos[n],year_indices_pos[n],spp_indices_pos[n]], obs_mat[spp_indices_pos[n],site_indices_pos[n]]);
}
