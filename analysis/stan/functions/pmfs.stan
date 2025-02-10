// Calculate the daily probability of reporting using parametric
// distributions up to the maximum observed delay.
// Adapted from https://github.com/epiforecasts/epinowcast
// (MIT License, copyright: epinowcast authors)

// taken from: https://github.com/epiforecasts/EpiNow2/blob/b3857a9accbe56222a9d816e7cf5d8473063cc32/inst/stan/functions/pmfs.stan

vector discretised_pmf(vector params, int n, int dist) {
  vector[n] lpmf;
  vector[n] upper_lcdf;
  if (dist == 0) {
    for (i in 1:n) {
      upper_lcdf[i] = lognormal_lcdf(i | params[1], params[2]);
    }
  } else if (dist == 1) {
    for (i in 1:n) {
      upper_lcdf[i] = gamma_lcdf(i | params[1], params[2]);
    }
  } else if (dist == 2) {
    for (i in 1:n) {
      upper_lcdf[i] = weibull_lcdf(i | params[1], params[2]);
    }
  } else {
    reject("Unknown distribution function provided.");
  }
  
  // discretise
  if (n > 1) {
    lpmf[1] = upper_lcdf[1];
    lpmf[2] = upper_lcdf[2];
    if (n > 2) {
      lpmf[3:n] = log_diff_exp(upper_lcdf[3:n], upper_lcdf[1:(n - 2)]);
    }
    // normalize
    lpmf = lpmf - log_sum_exp(upper_lcdf[(n - 1):n]);
  } else {
    lpmf[1] = 0;
  }
  return(exp(lpmf));
}

vector log_discretised_pmf(vector params, int n, int dist) {
  vector[n] lpmf;
  vector[n] upper_lcdf;
  if (dist == 0) {
    for (i in 1:n) {
      upper_lcdf[i] = lognormal_lcdf(i | params[1], params[2]);
    }
  } else if (dist == 1) {
    for (i in 1:n) {
      upper_lcdf[i] = gamma_lcdf(i | params[1], params[2]);
    }
  } else if (dist == 2) {
    for (i in 1:n) {
      upper_lcdf[i] = weibull_lcdf(i | params[1], params[2]);
    }
  } else {
    reject("Unknown distribution function provided.");
  }
  // discretise
  if (n > 1) {
    lpmf[1] = upper_lcdf[1];
    lpmf[2] = upper_lcdf[2];
    if (n > 2) {
      lpmf[3:n] = log_diff_exp(upper_lcdf[3:n], upper_lcdf[1:(n - 2)]);
    }
    // normalize
    lpmf = lpmf - log_sum_exp(upper_lcdf[(n - 1):n]);
  } else {
    lpmf[1] = 0;
  }
  return(lpmf);
}


// Function to compute the log of cdf of target distribution
real log_cdf_dist(vector params, real bound, int dist) {
  real log_cdf_val;
  
  if (dist == 0) {
    log_cdf_val = lognormal_lcdf(bound | params[1], params[2]);
  } else if (dist == 1) {
    log_cdf_val = gamma_lcdf(bound | params[1], params[2]);
  } else if (dist == 2) {
    log_cdf_val = weibull_lcdf(bound | params[1], params[2]);
  } else {
    reject("Unknown distribution function provided.");
  }
  
  return(log_cdf_val);
}

// Compute delay likelihood accounting for daily censoring (Park et al. 2024)
real daily_censoring_ll(real dt, vector pars, int dist) {
  real ll;
  
  if (dt <= 1) {
    ll = log_cdf_dist(pars, dt + 1, dist);
  } else {
    real a = log_cdf_dist(pars, dt + 1, dist);
    real b = log_cdf_dist(pars, dt - 1, dist);
    ll = log_diff_exp(a, b);  
  }
  
  return ll;
}

// Compute probability of exceedence
real exceedence_prob(real thresh, vector pars, int dist) {
  real prop_exceed = 1;
  
  if (dist == 0) {
    prop_exceed -= lognormal_cdf(thresh | pars[1], pars[2]);
  } else if (dist == 1) {
    prop_exceed -= gamma_cdf(thresh | pars[1], pars[2]);
  } else if (dist == 2) {
    prop_exceed -= weibull_cdf(thresh | pars[1], pars[2]);
  }
  
  return prop_exceed;
}

