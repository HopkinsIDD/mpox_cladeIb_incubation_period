//
// Inference of incubation period - Uvira Mpox
// This is the base model that only accounts for censoring

#include infer_incubation_period_v0_trunc2.stan
generated quantities {
  matrix[J, 3] prop_exceed;    // probatbility of exceedence
  
  {
    row_vector[3] thresh = [14, 21, 28];
    for (l in 1:J) {
      for (i in 1:3) {
        prop_exceed[l, i] = exceedence_prob(thresh[i], pars[, l], dist);
      }
    }
  }
}
