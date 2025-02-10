// 
// Inference of incubation period - Uvira Mpox
// This is the base model that accounts for censoring and truncation
// 

functions {
  #include functions/convolve.stan
  #include functions/pmfs.stan
  #include functions/delays.stan
  #include functions/utils.stan
}
data {
  // Numbers
  int<lower=1> N;    // number of observations (individuals x visits x possible contacts)
  int<lower=1> M;    // number of individuals x visits
  int<lower=1> J;    // number of groups for which the incubation period is inferred
  int<lower=0> L;    // number of covariates
  int<lower=0> O;
  int<lower=0> P;    // number of days of of simulated infections
  int<lower=0> Q;    // max number of days of reporting delay
  int<lower=0> N_delays;    // number of observed delays between onset and reporting 
  array[J] int<lower=0> N_observations;    // number of days of case observations
  array[J] int<lower=0> N_no_observations;    // number of days of case observations
  
  // Covariates
  matrix[M, L] X;    // covariate matrix
  
  // Data
  array[N] real<lower=0> t_s;       // time of symptom onset (scaled to days from study start)
  array[N] real<lower=0> t_e_r;     // upper bound for time of exposure (scaled to days from study start)
  array[N] real<lower=0> t_e_l;     // upper bound for time of exposure (scaled to days from study start)
  array[N] int<lower=0> t_s_int;     // upper bound for time of exposure (scaled to days from study start)
  array[N] int<lower=0,upper=1> censored;     // whether exposure window is known or not
  
  array[N] int<lower=1, upper=J> map_to_group;    // map of observations to the group
  array[M] int<lower=1, upper=J> map_patient_to_group;    // map of patients to the group
  array[M] int<lower=1, upper=N> starts;    // start of observations for each individual x visit
  array[M] int<lower=1, upper=N> ends;      // end of observations for each individual x visit
  array[M] int<lower=1> n_contacts;         // number of contacts for each invidual x visit
  
  array[P, J] int y;    // daily number of reported infections
  vector<lower=0>[N_delays] dt_reporting;    // daily number of reported infections
  array[N_delays] int<lower=0,upper=J> map_delays_group;
  array[J, max(N_observations)] int<lower=0,upper=P> ind_observations;         // number of contacts for each invidual x visit
  array[J, max(N_no_observations)] int<lower=0,upper=P> ind_no_observations;         // number of contacts for each invidual x visit
  
  // Prior on duration of lesions, assumed to follow a log-normal distribution
  // We here assume that the infectious period corresponds to the lesion periods
  real<lower=0> log_mu_lesions;
  real<lower=0> sigma_lesions;
  
  // Priors incubation period
  real<lower=0> sd_prior;
  real<lower=0> log_mu_prior;
  real<lower=0> sigma_sd_prior;
  real<lower=0> sigma_mu_prior;
  
  // Options for incubation period distribution
  int<lower=0, upper=2> dist;
  
  // Auxilliary data
  array[N] int<lower=0> n_steps_prior;    // number of days over which to integrate
  int<lower=1> n_days_prior;
  int<lower=1> n_q_gen;
  int<lower=0> n_target_q_gen;
  vector[n_target_q_gen] p_target;
  
  // Flags
  int<lower=0, upper=1> debug;
  int<lower=0, upper=1> use_data;
  int<lower=0, upper=1> truncate;
  int<lower=0, upper=1> incid_by_group;
  
  // Number of untruncated observations
  int<lower=0,upper=N> n_untrunc;
  
  // Auxilliary
  array[0] real x_r;
  array[0] int x_i;
  
}

transformed data {
  vector[N] log_prior_contacts;    // prior on who is the infecting contact  
  array[N] vector[max(n_steps_prior)] log_prior_exposures;    // prior on when exposure occurred
  int multigroup;
  int tot_delay = O+Q;
  int N_incid = P+tot_delay;
  array[P] int tot_y;    // sum of cases across groups
  int J_incid;
  
  if (J > 1) {
    multigroup = 1;
  } else {
    multigroup = 0;
  }
  
  if (J > 1 && incid_by_group == 1) {
    J_incid = J;
  } else {
    J_incid = 1;
  }
  
  for (i in 1:P) {
    tot_y[i] = sum(y[i, ]);
  }
  
  for (m in 1:M) {
    for (i in starts[m]:ends[m]) {
      // assuming uniform prior over contacts
      log_prior_contacts[i] = log(1./n_contacts[m]);
    }
  }
  
  for (i in 1:N) {
    int n = n_steps_prior[i];
    array[n] real ts;
    ts = linspaced_array(n, t_e_r[i] - n, t_e_r[i]);
    
    for (j in 1:n) {
      log_prior_exposures[i][j] = lognormal_lccdf(t_e_r[i] - ts[j] + 1| log_mu_lesions, sigma_lesions);
      if (is_nan(log_prior_exposures[i][j] )) {
        print("ll prior, i", i, " j", j, " ll", log_prior_exposures[i][j]);
      }
    }
    
    // Normalize
    log_prior_exposures[i][1:n] = log(softmax(log_prior_exposures[i][1:n]));
  }
  
  // Statements
  print("Running with flags: use_data:", use_data, "; dist:", dist);
}

parameters {
  vector[J] log_mu_std;                // location parameter of incubation period
  vector<lower=0>[J] sigma;    // scale paremeter of incubation period
  vector[multigroup] log_mu_mu;
  vector<lower=0,upper=5>[multigroup] sd_mu;
  real logit_phi;
  vector[J] mu_logit_lambda;
  matrix[P, J_incid] log_rho;      // latent infections in community (starting O days prior to first reported case)  
  vector[J_incid] mu_rho;
  vector<lower=0>[J_incid] sd_rho;
  vector[J] log_mu_reporting;
  vector<lower=0>[J] sigma_reporting;
  real<lower=0> inv_od;
}

transformed parameters {
  vector[J] log_mu;
  vector[J] mu;                // location parameter of incubation period
  vector[multigroup] mu_mu;                // location parameter of incubation period
  // map mean/sd to par1/par2 of distributions, which depend on each function
  matrix[2, J] pars;
  vector[use_data*M] log_lik;
  vector[J] log_cdf_trunc;    // account for truncation of reported contacts
  real phi = inv_logit(logit_phi);
  matrix[N_incid, J_incid] infections;
  matrix[O+1, J] disc_inc_period;
  matrix[O+1, J] log_disc_inc_period;
  matrix[Q+1, J] disc_reporting;
  matrix[tot_delay+1, J] disc_delay_pmf;
  matrix[P, J] reported_cases;
  vector[P] tot_reported_cases;
  vector[J] mu_reporting = exp(log_mu_reporting);
  
  // Compute mean of incubation periods
  if (multigroup == 1){
    log_mu = log_mu_mu[1] + log_mu_std * sd_mu[1];
    mu_mu = exp(log_mu_mu);
  } else {
    log_mu = log_mu_prior + log_mu_std * sd_prior;
  }
  
  mu = exp(log_mu);
  
  // Map from mean/sd to parameterization of target distribution
  for (j in 1:J) {
    pars[, j] = map_params(mu[j], sigma[j], dist, x_r, x_i);
  }
  
  // Make vector of infections
  for (j in 1:J_incid) {
    infections[, j] = append_row(rep_vector(0, tot_delay), exp(mu_rho[j] + log_rho[, j])) + 1e-3;
  }
  
  // Make disctretized pmfs and convolutions
  for (j in 1:J) {
    disc_inc_period[, j] = discretised_pmf(pars[, j], O+1, dist);
    log_disc_inc_period[, j] = log_discretised_pmf(pars[, j], O+1, dist);
    disc_reporting[, j] = discretised_pmf([log_mu_reporting[j], sigma_reporting[j]]', Q+1, 0);
    disc_delay_pmf[, j] = convolve_with_rev_pmf(disc_inc_period[, j], reverse(disc_reporting[, j]), tot_delay+1);
    
    // if (is_nan(disc_delay_pmf[1, j])) {
    //   print("nan delay, inc:", disc_inc_period[1, j], " report:", disc_reporting[1, j]);
    // }
    
    if (incid_by_group == 1) {
      reported_cases[, j] = convolve_to_report(infections[, j], reverse(disc_delay_pmf[, j]), tot_delay);
    } else {
      reported_cases[, j] = convolve_to_report(infections[, 1], reverse(disc_delay_pmf[, j]), tot_delay);
    }
  }
  
  // Compute total reported cases
  for (i in 1:P) {
    tot_reported_cases[i] = sum(reported_cases[i, ]);
    // if (is_nan(tot_reported_cases[i]) || tot_reported_cases[i] < 0) {
    //   print("nan cases i:", i, "inc:", disc_inc_period[1, 1], " report:", disc_reporting[1, 1], " lmureport:", log_mu_reporting[1], " sd:", sigma_reporting[1]);
    // }
  }
  
  // Compute the cdf at 21 days to account for truncation in questionnaire
  for (j in 1:J) {
    log_cdf_trunc[j] = log_cdf_dist(pars[, j], 21, dist);
  }
  
  if (use_data == 1) {
    for (m in 1:M) {
      int n = n_contacts[m];   // this particiant's number of reported contacts
      int cnt = 1;    // coutner for contacts for this participant
      vector[n] ll;   // placeholder for log-liks for each contact
      vector[O+1] ll_comm;   // placeholder for log-liks for exposures from community
      vector[2] ll_tot;
      real lambda = inv_logit(mu_logit_lambda[map_patient_to_group[m]]);
      
      // Integrate over all contacts
      for (i in starts[m]:ends[m]) {
        real ll_c;
        int l = map_to_group[i];    // index of incubation period group
        
        if (censored[i] == 1) {
          int k = n_steps_prior[i];          // Case of right censored incubation time
          vector[k] ll_exp;    // log_lik of exposures prior to most recent contact
          vector[k] ll_exp2 = rep_vector(-log(k), k);    // ll denominator to account for truncation 
          array[k] real ts;    // times of possible exposures to integrate over
          
          ts = linspaced_array(k, t_e_r[i] - k, t_e_r[i]);
          
          for (j in 1:k) {
            real dt = t_s[i] - ts[j];
            real bound = t_e_r[i] - ts[j] + 21;
            real this_ll;
            
            // Compute likelihood accounting for daily censoring 
            this_ll = daily_censoring_ll(dt, pars[, l], dist);
            
            ll_exp[j] = log_prior_exposures[i][j] + this_ll;
            
            // Denominator for truncation
            if ((t_s[i] - t_e_r[i]) <= 21) {
              if (bound < 50) {
                ll_exp2[j] = log_prior_exposures[i][j] + log_cdf_dist(pars[, l], bound, dist);
              } else {
                ll_exp2[j] = log_prior_exposures[i][j];
              }
            }
          }
          
          // Likelihood for this contact
          ll_c = log_sum_exp(ll_exp);
          
          // Account for truncation 
          if ((t_s[i] - t_e_r[i]) <= 21 && truncate == 1) {
            ll_c = log_mix(phi, ll_c, ll_c - log_sum_exp(ll_exp2));
          }
          
        } else {
          // Case where single exposure to contact, only account for daily censoring
          real dt = t_s[i] - t_e_l[i];
          
          // Compute likelihood accounting for daily censoring
          ll_c  = daily_censoring_ll(dt, pars[, l], dist);
          
          // Right-truncation for reporting contacts in the past 21 days
          if (dt <= 21 && truncate == 1) {
            ll_c = log_mix(phi, ll_c, ll_c - log_cdf_trunc[l]);
          }
          
          // if (is_nan(ll_c)) {
          //   print("ll_c nan at dt:", dt, " pars1:", pars[1, l], " pars2:",  pars[2, l],  " ll", log_cdf_trunc[l]);
          // }
        }
        
        // this contact's log-lik
        ll[cnt] = log_prior_contacts[i] + ll_c;
        cnt += 1;
      }
      
      // Integrate exposures from the community
      {
        int j = starts[m];
        int l = map_to_group[j];    // index of incubation period group
        int tL = t_s_int[j] + tot_delay - O;
        int tR = t_s_int[j] + tot_delay;
        vector[O+1] eta;
        
        if (incid_by_group == 1) {
          eta = infections[tL:tR, l]/sum(infections[tL:tR, l] + 1e-3);
        } else {
          eta = infections[tL:tR, 1]/sum(infections[tL:tR, 1] + 1e-3);
        }
        
        for (i in 1:(O+1))  {
          int k = O - i + 2;
          ll_comm[i] = log(eta[i]) + log_disc_inc_period[k, l];
        }
      }
      
      // Sum over contacts and community
      log_lik[m] = log_mix(lambda, log_sum_exp(ll), log_sum_exp(ll_comm));
    }
  }
}

model {
  
  if (use_data == 1) {
    target += sum(log_lik);
  }
  
  // Priors
  if (J == 1) {
    log_mu ~ normal(log_mu_prior, sd_prior);
  } else {
    log_mu_mu ~ normal(log_mu_prior, sd_prior);
    log_mu_std ~ std_normal();
    sd_mu ~ normal(0, .5);
  }
  sigma ~ normal(sigma_mu_prior, sigma_sd_prior);
  logit_phi ~ normal(-1, 1);
  target += binomial_lccdf(n_untrunc | N, phi);
  
  if (J > 1) {
    // Observed cases
    for (j in 1:J) {
      array[N_observations[j]] int ind_obs = ind_observations[j, 1:N_observations[j]];
      array[N_no_observations[j]] int ind_no_obs = ind_no_observations[j, 1:N_no_observations[j]];
      
      y[ind_obs, j] ~ neg_binomial_2(reported_cases[ind_obs, j], 1/inv_od);
      0 ~ poisson(reported_cases[ind_no_obs, j]);
    }
  } else {
    tot_y[ind_observations[1]] ~ neg_binomial_2(tot_reported_cases[ind_observations[1]], 1/inv_od);
    
    0 ~ poisson(tot_reported_cases[ind_no_observations[1]]);
  }
  
  for (j in 1:J_incid) {
    target += autocorrPrior(P, 1/sd_rho[j]^2, log_rho[, j], 1);
  }
  
  // Dispersion of reporting
  inv_od ~ normal(0, .5);
  
  // Infections
  sd_rho ~ normal(0, .25);
  for (i in 1:J_incid) {
    log_rho[1, i] ~ normal(-2, 1);
  }
  mu_rho ~ normal(0, 1);
  
  // Probability of infection in contacts
  mu_logit_lambda ~ normal(2, .4);
  
  // Reporting delays
  dt_reporting ~ lognormal(log_mu_reporting[map_delays_group], sigma_reporting[map_delays_group]);
  log_mu_reporting ~ normal(log(5), .3);
  sigma_reporting ~ normal(0, .5);
  
}

