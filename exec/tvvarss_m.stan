functions {
  \\ mapB: maps b onto [lo,hi]
  \\ inv_logit(b) = mapB(b, 0, 1)
  real mapB(real b, real lo, real up) {
    (up - lo) / (1 + exp(-b)) + lo
  }
}
data {
  // scalars
  int<lower=0> n_year;                      // # years
  int<lower=0> n_site;                      // # sites
  int<lower=0> n_spp;                       // # species
  int<lower=0> n_process;                   // # of total processes
  int<lower=0> n_q;                         // # proc variances = max(shared_q)
  int<lower=0> n_r;                         // # obs variances = max(shared_r)
  int<lower=0> n_u;                         // # trends = max(shared_u)
  int<lower=0> n_pos;                       // # non-NA values in the data
  int<lower=0> fit_dynamicB;                // indicator of fit dynamic B: Y/N = 1/0
  int<lower=0> est_trend;                   // indicator of bias in RW: Y/N = 1/0
  int<lower=0> demean;                      // indicator of demean proc: Y/N = 1/0
  int family;                               // indicator of distn for obs data
  // vectors
  int<lower=0> process[n_site+1];           // map of sites:processes
  int<lower=0> b_diag[n_spp*n_spp];         // indicators of on-diagonal: Y/N = 2/1
  int<lower=0> b_indices[(n_spp*n_spp)];    // indices for B values
  int<lower=0> row_indices[(n_spp*n_spp)];  // row indices for B values
  int<lower=0> col_indices[(n_spp*n_spp)];  // col indices for B values
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
parameters {
  vector[(n_spp*n_spp)] vecBdev[n_year]; // elements accessed [n_year,n_spp]
  real<lower=0> sigma_rw_pars[2]; // sds for random walk
  matrix[n_year,n_spp] x[n_process]; // unobserved states
  real<lower=0> resid_process_sd[n_q]; // residual sds
  real<lower=0> obs_sd[n_r]; // residual sds
  real u[n_u]; // trends
}
transformed parameters {
  vector<lower=0>[(n_spp*n_spp)] sigma_rw;
  matrix<lower=0>[n_spp, n_process] resid_process_mat;
  matrix<lower=0>[n_spp, n_site] obs_mat;
  vector<lower=-20,upper=20>[(n_spp*n_spp)] vecB[n_year];
  matrix[n_spp, n_process] u_mat;
  matrix[n_spp,n_spp] B[(n_year-1)]; // B matrix, accessed as n_year, n_spp, n_spp
  matrix[n_year,n_spp] pred[n_process]; // predicted unobserved states

  for(i in 1:(n_spp*n_spp)) {
    sigma_rw[i] = sigma_rw_pars[b_diag[i]];
  }
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
  // for(i in 1:(n_spp*n_spp)) {
    // vecB[1,i] = vecBdev[1,i]; // first time step, prior in model {}
  // }
  vecB[1] = vecBdev[1];
  for(t in 2:n_year) {
    // fill in B matrix, shared across sites

    // estimated interactions
    for(i in 1:(n_spp*n_spp)) {
      // top down (td) interactions
      B[t-1,row_indices[i],col_indices[i]] = (b_limits[b_indices[i],2] - b_limits[b_indices[i],1]) * vecB[t-1,i] + b_limits[b_indices[i],1];
      // if user wants constant b matrix, vecB[t] = vecB[t-1]
      vecB[t,i] = vecB[t-1,i] + (fit_dynamicB)*vecBdev[t-1,i]; // random walk in b elements
    }

    //for(i in 1:(n_spp*n_spp)) {
    //  B[t-1,row_indices[i],col_indices[i]] = vecB[t-1,i];
    //}

    // do projection to calculate predicted values. modify code depending on whether
    // predictions should be demeaned before projected, and whether or not trend included.
    for(s in 1:n_process) {
     if(est_trend == 0) {
      if(demean==0) {pred[s,t,] = x[s,t-1,] * B[t-1,,];}
      if(demean==1) {pred[s,t,] = (x[s,t-1,]- u_mat[,s]') * B[t-1,,];}
     }
     if(est_trend == 1) {
      if(demean==0) {pred[s,t,] = x[s,t-1,] * B[t-1,,] + u_mat[,s]';}
      if(demean==1) {pred[s,t,] = (x[s,t-1,]-u_mat[,s]') * B[t-1,,] + u_mat[,s]';}
     }
    }
  }
}
model {
  sigma_rw_pars[1] ~ student_t(5,0,1);
  sigma_rw_pars[2] ~ student_t(5,0,1);
  //vecB[1] ~ normal(0, 3); // prior for first state
  for(t in 1:n_year) {
    //vecB[t] ~ normal(vecB[t-1], sigma_rw); // vectorized random in B
    vecBdev[t] ~ normal(0, 10); // vectorized random in B // removed (0, sigma_rw)
  }
  // prior on first time step
  for(site in 1:n_process) {
    for(spp in 1:n_spp) {
      x[site,1,spp] ~ normal(x0[site,spp],1);
    }
  }
  // process model for remaining sites
  for(t in 2:n_year) {
   for(site in 1:n_process) {
    for(spp in 1:n_spp) {
      x[site,t,spp] ~ normal(pred[site,t,spp], resid_process_mat[spp,site]);
    }
   }
  }

  for(i in 1:n_q) {
    // prior on process standard deviations
    resid_process_sd[i] ~ student_t(5,0,1);
  }
  for(i in 1:n_r) {
    obs_sd[i] ~ student_t(5,0,1);
  }
  for(i in 1:n_u) {
    // prior on trends
    u[i] ~ normal(0,1);
  }

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
  vector[n_pos] log_lik;
  // for use in loo() package
  if(family==1) for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | x[site_indices_pos[n],year_indices_pos[n],spp_indices_pos[n]], obs_mat[spp_indices_pos[n],site_indices_pos[n]]);
  if(family==2) for (n in 1:n_pos) log_lik[n] = bernoulli_lpmf(y_int[n] | inv_logit(x[site_indices_pos[n],year_indices_pos[n],spp_indices_pos[n]]) );
  if(family==3) for (n in 1:n_pos) log_lik[n] = poisson_lpmf(y_int[n] | exp(x[site_indices_pos[n],year_indices_pos[n],spp_indices_pos[n]]) );
  if(family==4) for (n in 1:n_pos) log_lik[n] = gamma_lpdf(y[n] | obs_mat[spp_indices_pos[n],site_indices_pos[n]], obs_mat[spp_indices_pos[n],site_indices_pos[n]] ./ x[site_indices_pos[n],year_indices_pos[n],spp_indices_pos[n]]);
  if(family==5) for (n in 1:n_pos) log_lik[n] = lognormal_lpdf(y[n] | x[site_indices_pos[n],year_indices_pos[n],spp_indices_pos[n]], obs_mat[spp_indices_pos[n],site_indices_pos[n]]);
}
