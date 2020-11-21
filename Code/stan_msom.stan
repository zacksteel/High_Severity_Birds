// Full Multi-species occupancy model

functions {
  
  real lp_observed(int y, int J, real logit_psi, real logit_p) {
    return log_inv_logit(logit_psi)
      + binomial_logit_lpmf(y | J, logit_p);

  }

  real lp_unobserved(int J, real logit_psi, real logit_p) {
    return log_sum_exp(lp_observed(0, J, logit_psi, logit_p),
      log1m_inv_logit(logit_psi));
  }
}
                          

data {
  int<lower=1> S;                 // Number of species
  int<lower=1> W;                 // Number of wildfires
  int<lower=1> L;                 // Number of locations (sites)
  int<lower=1> P;                 // Number of Patches
  int<lower=1> M;                 // Number of surveys (sites x years)
  int<lower=1> J;                 // Number of temporal replications

  int<lower=0,upper=1> y[M, J, S];  // Observation at site M, visit J, for species S
  int locs[M];                    // site/location index
  int patches[M];                 // patch index
  int fires[M];                   // Fire index
  
  vector[M] dist;                 // Distance Covariate
  vector[M] ysf;                  // Year since fire covariate
  vector[M] ysfsq;                // YSF quadratic
  vector[M] dy;                   // dist x ysf interaction term
  vector[M] parea;                // patch area
  vector[M] elev;                 // elevation
  vector[M] elevsq;               // elevation quadratic
  vector[M] lat;                  // latitude
  vector[M] snagba;               // snag basal area
  vector[M] shrub;                // shrub cover
  matrix[M, J] jday;              // julian day 
  matrix[M, J] jdaysq;            // jday quadratic
  
  int last[M];                    // indicator of how many visits were made per site
}

parameters {
  
  //Intercepts
  
  real beta_mu; // hyperparameter mean
  real<lower=0> beta_sd; // hyperparameter sd
  vector[S] beta_tilde; // non-centered presence coefficient

  real alpha_mu; 
  real<lower=0> alpha_sd; 
  vector[S] alpha_tilde;
  
  // Species-specific random intercepts
  
  // // Location (site)
  vector<lower=0>[S] sigma_b_loc;
  matrix[S,L] tilde_b_loc;
  
  // Patch
  vector<lower=0>[S] sigma_b_patch;
  matrix[S,P] tilde_b_patch;
  
  // Fire
  vector<lower=0>[S] sigma_b_fire;
  matrix[S,W] tilde_b_fire;

  
  //Slopes
  real mu_b_dist;   // distance; community mean
  real<lower=0> sigma_b_dist; //community sigma
  vector[S] tilde_b_dist;

  real mu_b_ysf; //ysf
  real<lower=0> sigma_b_ysf;
  vector[S] tilde_b_ysf;

  real mu_b_ysfsq; //ysf squared
  real<lower=0> sigma_b_ysfsq;
  vector[S] tilde_b_ysfsq;

  real mu_b_dy; //distance * ysf interaction
  real<lower=0> sigma_b_dy;
  vector[S] tilde_b_dy;
  
  real mu_b_parea; //high severity area residual
  real<lower=0> sigma_b_parea;
  vector[S] tilde_b_parea;

  real mu_b_elev; //elevation
  real<lower=0> sigma_b_elev;
  vector[S] tilde_b_elev;

  real mu_b_elevsq; //elevation squared
  real<lower=0> sigma_b_elevsq;
  vector[S] tilde_b_elevsq;
  
  real mu_b_lat; //latitude
  real<lower=0> sigma_b_lat;
  vector[S] tilde_b_lat;

  real mu_a_snagba;
  real<lower=0> sigma_a_snagba;
  vector[S] tilde_a_snagba;

  real mu_a_shrub;
  real<lower=0> sigma_a_shrub;
  vector[S] tilde_a_shrub;

  real mu_a_jday; // julian day
  real<lower=0> sigma_a_jday;
  vector[S] tilde_a_jday;

  real mu_a_jdaysq; // jday squared
  real<lower=0> sigma_a_jdaysq;
  vector[S] tilde_a_jdaysq;
  
}

transformed parameters {
  
  // non-centered global intercepts
  vector[S] beta = beta_mu + beta_sd * beta_tilde;
  vector[S] alpha = alpha_mu + alpha_sd * alpha_tilde;
  
  // non-centered slopes
  vector[S] b_dist = mu_b_dist + sigma_b_dist * tilde_b_dist;
  vector[S] b_ysf = mu_b_ysf + sigma_b_ysf * tilde_b_ysf;
  vector[S] b_ysfsq = mu_b_ysfsq + sigma_b_ysfsq * tilde_b_ysfsq;
  vector[S] b_dy = mu_b_dy + sigma_b_dy * tilde_b_dy;
  vector[S] b_parea = mu_b_parea + sigma_b_parea * tilde_b_parea;
  vector[S] b_elev = mu_b_elev + sigma_b_elev * tilde_b_elev;
  vector[S] b_elevsq = mu_b_elevsq + sigma_b_elevsq * tilde_b_elevsq;
  vector[S] b_lat = mu_b_lat + sigma_b_lat * tilde_b_lat;

  vector[S] a_snagba = mu_a_snagba + sigma_a_snagba * tilde_a_snagba;
  vector[S] a_shrub = mu_a_shrub + sigma_a_shrub * tilde_a_shrub;
  vector[S] a_jday = mu_a_jday + sigma_a_jday * tilde_a_jday;
  vector[S] a_jdaysq = mu_a_jdaysq + sigma_a_jdaysq * tilde_a_jdaysq;
  
  
  // species-specific random intercepts
  matrix[S,W] b_fire;
  matrix[S,P] b_patch;
  matrix[S,L] b_loc;
  
  // Below Moved from model block
  
  real logit_psi[M,S];          // Logit occupancy probability, each species
  real logit_p[M,S];            // Logit detection probability, each species

  for (s in 1:S) {
    b_loc[s] = sigma_b_loc[s] * tilde_b_loc[s]; // mu set to zero
    b_patch[s] = sigma_b_patch[s] * tilde_b_patch[s];
    b_fire[s] = sigma_b_fire[s] * tilde_b_fire[s];
  }
  
  for (s in 1:S) {              // Loop through each species
    for (i in 1:M) {            // Loop through each site
      logit_psi[i,s] = beta[s] + 
      b_loc[s,locs[i]] + b_patch[s,patches[i]] + b_fire[s,fires[i]] + //species-level intercepts
      b_dist[s] * dist[i] + b_ysf[s] * ysf[i] + b_ysfsq[s] * ysfsq[i] + 
      b_dy[s] * dy[i] + 
      b_parea[s] * parea[i] +
      b_elev[s] * elev[i] + b_elevsq[s] * elevsq[i] + b_lat[s] * lat[i];
    }
  }

  // linear model for detection
  for (s in 1:S) {
    for (i in 1:M) {
      for (j in 1:J) {
        logit_p[i,s] = alpha[s] + 
        a_snagba[s] * snagba[i] + a_shrub[s] * shrub[i] +
        a_jday[s] * jday[i,j] + a_jdaysq[s] * jdaysq[i,j];
      }
    }
  }

}

model {

  // Priors
  // Intercepts
  beta_mu ~ normal(0, 1); // implies beta ~ normal(beta_mu, beta_sd)
  beta_sd ~ normal(0, 1);
  beta_tilde ~ normal(0, 1); 
	
  alpha_mu ~ normal(0, 1);
  alpha_sd ~ normal(0, 1);
  alpha_tilde ~ normal(0, 1);
  
  for (s in 1:S) {
    sigma_b_loc[s] ~ normal(0, 1);
    sigma_b_patch[s] ~ normal(0, 1);
    sigma_b_fire[s] ~ normal(0, 1);
    tilde_b_loc[s] ~ normal(0, 1);
    tilde_b_patch[s] ~ normal(0, 1);
    tilde_b_fire[s] ~ normal(0, 1);
  }
  
  tilde_b_dist ~ normal(0, 1);
  tilde_b_ysf ~ normal(0, 1);
  tilde_b_ysfsq ~ normal(0, 1);
  tilde_b_dy ~ normal(0, 1);
  tilde_b_parea ~ normal(0, 1);
  tilde_b_elev ~ normal(0, 1);
  tilde_b_elevsq ~ normal(0, 1);
  tilde_b_lat ~ normal(0, 1);

  tilde_a_snagba ~ normal(0, 1);
  tilde_a_shrub ~ normal(0, 1);
  tilde_a_jday ~ normal(0, 1);
  tilde_a_jdaysq ~ normal(0, 1);

  sigma_b_dist ~ normal(0, 1);
  sigma_b_ysf ~ normal(0, 1);
  sigma_b_ysfsq ~ normal(0, 1);
  sigma_b_dy ~ normal(0, 1);
  sigma_b_parea ~ normal(0, 1);
  sigma_b_elev ~ normal(0, 1);
  sigma_b_elevsq ~ normal(0, 1);
  sigma_b_lat ~ normal(0, 1);
 
  sigma_a_snagba ~ normal(0, 1);
  sigma_a_shrub ~ normal(0, 1);
  sigma_a_jday ~ normal(0, 1);
  sigma_a_jdaysq ~ normal(0, 1);

  mu_b_dist ~ normal(0, 1);
  mu_b_ysf ~ normal(0, 1);
  mu_b_ysfsq ~ normal(0, 1);
  mu_b_dy ~ normal(0, 1);
  mu_b_parea ~ normal(0, 1);
  mu_b_elev ~ normal(0, 1);
  mu_b_elevsq ~ normal(0, 1);
  mu_b_lat ~ normal(0, 1);

  mu_a_snagba ~ normal(0, 1);
  mu_a_shrub ~ normal(0, 1);
  mu_a_jday ~ normal(0, 1);
  mu_a_jdaysq ~ normal(0, 1);

  // Likelihood
  for (s in 1:S) {
    for (i in 1:M) {
      if (sum(y[i,,s]) > 0) 
        target += lp_observed(sum(y[i,,s]), last[i], logit_psi[i,s], logit_p[i,s]);
      else
        target += lp_unobserved(last[i], logit_psi[i,s], logit_p[i,s]);        
    }
  }
}

// generated quantities {
// 
//   real psi_con[M,S];  // prob occupied not conditional on data
//   real z[M,S];         // occupancy indicator, 0/1
// 
//   for (s in 1:S) {
//     for (i in 1:M) {
//       if (sum(y[i,,s]) > 0) {
//         psi_con[i,s] = lp_observed(sum(y[i,,s]), last[i], logit_psi[i,s], logit_p[i,s]);
//         z[i,s] = 1;
//       }
//       else {
//         psi_con[i,s] = lp_unobserved(last[i], logit_psi[i,s], logit_p[i,s]);
//         z[i,s] = bernoulli_logit_rng(psi_con[i,s]);
//       }
// 
//     }
//   }
// }

