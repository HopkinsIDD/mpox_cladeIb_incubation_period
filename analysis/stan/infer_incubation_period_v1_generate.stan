//
// Inference of incubation period - Uvira Mpox
// This is the base model that only accounts for censoring

#include infer_incubation_period_v1_trunc2.stan
generated quantities {
  array[P, J] int gen_reported_cases;
  array[P] int gen_tot_reported_cases;
  matrix[J, 3] prop_exceed;    // probatbility of exceedence
  
  {
    row_vector[3] thresh = [14, 21, 28];
    for (l in 1:J) {
      for (i in 1:3) {
        prop_exceed[l, i] = exceedence_prob(thresh[i], pars[, l], dist);
      }
    }
  }
  
  for (i in 1:P) {
    for (j in 1:J) {
      gen_reported_cases[i, j] = neg_binomial_2_rng(reported_cases[i, j], 1/inv_od);
    }
    gen_tot_reported_cases[i] = neg_binomial_2_rng(tot_reported_cases[i], 1/inv_od);
  }
}
