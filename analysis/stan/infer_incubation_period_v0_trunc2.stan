// 
// Inference of incubation period - Uvira Mpox
// This is the base model that accounts for censoring and truncation
// 

functions {
  #include functions/utils.stan
  #include functions/pmfs.stan
}
data {
  // Numbers
  int<lower=1> N;    // number of observations (individuals x visits x possible contacts)
  int<lower=1> M;    // number of individuals x visits
  int<lower=1> J;    // number of groups for which the incubation period is inferred
  int<lower=0> L;    // number of covariates
  int<lower=0> O;
  // Covariates
  matrix[M, L] X;    // covariate matrix
  
  // Data
  array[N] real<lower=0> t_s;       // time of symptom onset (scaled to days from study start)
  array[N] real<lower=0> t_e_r;     // upper bound for time of exposure (scaled to days from study start)
  array[N] real<lower=0> t_e_l;     // upper bound for time of exposure (scaled to days from study start)
  array[N] int<lower=0,upper=1> censored;     // whether exposure window is known or not
  
  array[N] int<lower=1, upper=J> map_to_group;    // map of observations to the group
  array[M] int<lower=1, upper=N> starts;    // start of observations for each individual x visit
  array[M] int<lower=1, upper=N> ends;      // end of observations for each individual x visit
  array[M] int<lower=1> n_contacts;         // number of contacts for each invidual x visit
  
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
  int truncate;
  
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
  
  if (J > 1) {
    multigroup = 1;
  } else {
    multigroup = 0;
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
    
    if (n > 0) {
      for (j in 1:n) {
        log_prior_exposures[i][j] = lognormal_lccdf(t_e_r[i] - ts[j] + 1| log_mu_lesions, sigma_lesions);
        if (is_nan(log_prior_exposures[i][j] )) {
          print("ll prior, i", i, " j", j, " ll", log_prior_exposures[i][j]);
        }
      }
      
      // Normalize
      log_prior_exposures[i][1:n] = log(softmax(log_prior_exposures[i][1:n]));
      
      if (i == 1) {
        print(log_prior_exposures[i]);
      }
    }
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
  
  // Compute the cdf at 21 days to account for truncation in questionnaire
  for (j in 1:J) {
    log_cdf_trunc[j] = log_cdf_dist(pars[, j], 21, dist);
  }
  
  if (use_data == 1) {
    for (m in 1:M) {
      int n = n_contacts[m];   // this particiant's number of reported contacts
      int cnt = 1;    // coutner for contacts for this participant
      vector[n] ll;   // placeholder for log-liks for each contact
      
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
          
          if (is_nan(ll_c)) {
            print("ll_c nan at dt:", dt, " pars1:", pars[1, l], " pars2:",  pars[2, l],  " ll", log_cdf_trunc[l]);
          }
        }
        
        // this contact's log-lik
        ll[cnt] = log_prior_contacts[i] + ll_c;
        cnt += 1;
      }
      
      // Sum over contacts
      log_lik[m] = log_sum_exp(ll);
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
}

