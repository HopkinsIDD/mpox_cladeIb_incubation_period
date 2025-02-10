// Utility functions for stan code

// Solver for weibull parameters
// from https://math.stackexchange.com/questions/1769765/weibull-distribution-from-mean-and-variance-to-shape-and-scale-factor
//  @param y vector of shape of weibull distribution
//  @param theta vector of mean and sd
//  @param x_r 
//  @param x_i
vector weibull_system(vector y, vector theta, data array[] real x_r, array[] int x_i) {
  vector[1] z;     
  real mu = theta[1];    // mean
  real sd2 = theta[2]^2;    // variance
  real alpha = exp(y[1]);    // shape parameter of weibull in stan function, k in link
  
  z[1] = sd2/mu^2 - tgamma(1+2/alpha)/tgamma(1+1/alpha)^2 + 1;
  return z;
}


// Convert mean/sd parametrization of weibull distribution to shape and scale
vector map_weibull_params(real mu, real sigma, data array[] real x_r, array[] int x_i) {
  vector[2] pars;
  vector[1] guess = [0]';
  vector[1] alpha;
  
  // checks
  alpha = exp(solve_powell(weibull_system, guess, [mu, sigma]', x_r, x_i));
  pars[1] = alpha[1];
  pars[2] = mu/(tgamma(1 + 1/alpha[1]));
  
  // print("pars: alpha ", pars[1], " sigma ", pars[2]);
  return pars;
}

// Convert mean/sd parametrization of gamma distribution to location and scale
vector map_gamma_params(real mu, real sigma) {
  vector[2] pars;
  
  // Location (alpha)
  pars[1] = (mu/sigma)^2;
  // Scale (beta)
  pars[2] = mu/sigma^2;
  
  return pars;
}

// Convert mean/sd parametrization of lognormal distribution to location and scale
vector map_lognormal_params(real mu, real sigma) {
  vector[2] pars;
  
  // Location log of the mean
  pars[1] = log(mu);
  // Scale sd on the log scale
  pars[2] = sigma;
  
  return pars;
}

// Returns mapped parameters
vector map_params(real mu, real sigma, int dist, data array[] real x_r, array[] int x_i) {
  vector[2] pars;
  
  if (dist == 0) {
    pars = map_lognormal_params(mu, sigma);
  } else if (dist == 1) {
    pars = map_gamma_params(mu, sigma);
  } else if (dist == 2) {
    pars = map_weibull_params(mu, sigma, x_r, x_i);
  }
  
  return pars;
}

// Prior for random walk
real autocorrPrior(int n, real tau, vector x, int order) {
  real ll;
  
  if (order == 1) {
    // First order autoregressive model
    ll = (n-1.0)/2.0 * log(tau) - tau/2 * (dot_self(x[2:n] - x[1:(n-1)]));
    // Soft sum to 0 constraint
    ll += normal_lpdf(sum(x) | 0, 0.001 * n);
  } else {
    vector[n] ix;
    for (i in 1:n) {
      ix[i] = i * x[i];
    }
    // Second order autoregressive model
    ll = (n-2.0)/2.0 * log(tau) - tau/2 * (dot_self(x[1:(n-2)] - 2*x[2:(n-1)] + x[3:n]));
    // Soft sum to 0 constraint
    ll += normal_lpdf(sum(x) | 0, 0.001 * n);
    // Soft sum to 0 constraint on sum_i ix
    ll += normal_lpdf(sum(ix) | 0, 0.001 * n);
  }
  
  return(ll);
}

