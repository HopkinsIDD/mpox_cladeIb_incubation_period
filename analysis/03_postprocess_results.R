# This script makes the figures and table for the incubation period analysis


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(posterior)
library(cowplot)
library(loo)
library(optparse)
library(here)
library(furrr)

source("analysis/utils.R")

# User-supplied options (if any)
option_list <- list(
  make_option(c("-r", "--redo_processing"), 
              default = TRUE, action ="store", type = "logical", 
              help = "Redo data processing.")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# Session for furrr
plan(multisession, workers = 10)

# Files
param_estimates_file <- here("generated_data/all_param_estimates.rds")
cdfs_file <- here("generated_data/all_cdfs.rds")
quantile_traj_file <- here("generated_data/all_quantile_traj.rds")
target_quantiles_file <- here("generated_data/all_target_quantiles.rds")
loo_file <- here("generated_data/all_loo.rds")
prop_exceed_file <- here("generated_data/all_prop_exceed.rds")
prop_exceed_draws_file <- here("generated_data/all_prop_exceed_draws.rds")

target_quantile_draws_file <- here("generated_data_quantile_draws.rds")

res_dir <- here("generated_data/idm_fits_Jan_30/")

# Load the results --------------------------------------------------------

# Stan fits
fit_files <- dir(res_dir, pattern = "fit", full.names = T) %>% 
  str_subset("dist-") %>% 
  str_subset("usedata-TRUE") %>% 
  str_subset("dic", negate = T)

overall_fits <- fit_files %>% 
  str_subset("group-none")

compare_fits <- fit_files %>% 
  str_subset("group-none")

# Stan fits
genquant_files <- dir(res_dir, pattern = "genquant", full.names = T) %>% 
  str_subset("dist-") %>% 
  str_subset("usedata-TRUE") %>% 
  str_subset("dic", negate = T)

overall_genquant_files <- genquant_files %>% 
  str_subset("group-none")

this_globals <- c("custom_quantile2", "custom_summaries", "cri_interval", "parse_inc_period_stan_filename")

# Mean/sd estimates -------------------------------------------------------
if (!file.exists(param_estimates_file) | opt$redo_processing) {
  cat("---- Extracting mean and sigma \n")
  
  param_estimates <- map_postprocess(
    files = fit_files,
    fun = postprocess_mu_sigma,
    .options = furrr_options(globals = c(this_globals, "postprocess_mu_sigma"),
                             packages = c("stringr", "dplyr"))
  )
  
  saveRDS(param_estimates, file = param_estimates_file)
} else {
  param_estimates <- readRDS(param_estimates_file)
}

# CDFS --------------------------------------------------------------------
if (!file.exists(cdfs_file) | opt$redo_processing) {
  cat("---- Extracting CDFs \n")
  
  inc_period_cdf <- map_postprocess(files = overall_fits,
                                    fun = postprocess_cdfs)
  
  saveRDS(inc_period_cdf, file = cdfs_file)
  
} else {
  inc_period_cdf <- readRDS(cdfs_file)
}

# Individual CDF draws --------------------------------------------------------------------
if (!file.exists(quantile_traj_file) | opt$redo_processing) {
  cat("---- Extracting quantile traj \n")
  
  inc_period_q_traj <- map_postprocess(files = overall_fits,
                                       fun = postprocess_quantiles_traj)
  
  saveRDS(inc_period_q_traj, file = quantile_traj_file)
  
} else {
  inc_period_q_traj <- readRDS(quantile_traj_file)
}

# Target quantiles --------------------------------------------------------
if (!file.exists(target_quantiles_file) | opt$redo_processing) {
  cat("---- Extracting quantiles \n")
  
  inc_period_q <- map_postprocess(files = fit_files,
                                  fun = postprocess_quantiles,
                                  .options = furrr_options(globals = c(this_globals, "postprocess_quantiles"),
                                                           packages = c("stringr", "dplyr")))
  
  saveRDS(inc_period_q, file = target_quantiles_file)
  
} else {
  inc_period_q <- readRDS(target_quantiles_file)
}

# Target quantile draws for stacking --------------------------------------
if (!file.exists(target_quantile_draws_file) | opt$redo_processing) {
  cat("---- Extracting quantile draws \n")
  
  inc_period_q_draws <- map_postprocess(files = overall_fits,
                                        fun = postprocess_quantiles_target_draws)
  
  saveRDS(inc_period_q_draws, file = target_quantile_draws_file)
  
} else {
  inc_period_q_draws <- readRDS(target_quantile_draws_file)
}

# Exceedence probs --------------------------------------------------------------------
if (!file.exists(prop_exceed_file) | opt$redo_processing) {
  cat("---- Exceedence probabilities \n")
  
  prop_exceed <- map_postprocess(files = overall_genquant_files,
                                 fun = postprocess_prop_exceed,
                                 .options = furrr_options(globals = c(this_globals, "postprocess_prop_exceed"),
                                                          packages = c("stringr", "dplyr")))
  
  saveRDS(prop_exceed, file = prop_exceed_file)
  
} else {
  prop_exceed <- readRDS(prop_exceed_file)
}

# Exceedence probs draws --------------------------------------------------------------------
if (!file.exists(prop_exceed_draws_file) | opt$redo_processing) {
  cat("---- Exceedence probabilities \n")
  
  prop_exceed_draws <- map_postprocess(files = overall_genquant_files,
                                       fun = postprocess_prop_exceed_draws)
  
  saveRDS(prop_exceed_draws, file = prop_exceed_draws_file)
  
} else {
  prop_exceed <- readRDS(prop_exceed_draws_file)
}

# Model comparison --------------------------------------------------------
if (!file.exists(loo_file) | opt$redo_processing) {
  cat("---- Extracting loo \n")
  
  loo_res <- map_postprocess(files = fit_files,
                             fun = postprocess_loo)
  
  saveRDS(loo_res, file = loo_file)
  
} else {
  loo_res <- readRDS(loo_file)
}

loo_compare <- loo_res %>% 
  filter(grouping == "none") %>%
  group_by(test_subset, symp, grouping) %>% 
  arrange(model_version, distribution) %>% 
  group_modify(function(x, y) {
    if (nrow(x) > 1) {
      
      print(y)
      res <- loo_compare(x$loo[1:nrow(x)])
      
      # Compute stacking weights
      w <- loo_model_weights(x$loo[1:nrow(x)], method = "stacking") %>% 
        as.data.frame() %>% 
        as_tibble()
      
      res %>% 
        as_tibble() %>% 
        mutate(across(everything(), as.numeric),
               id = str_extract(row.names(res), "[0-9]") %>% as.numeric()) %>% 
        inner_join(x %>% mutate(id = row_number(), weight = w$x)) %>% 
        mutate(rank = row_number(),
               significant = abs(elpd_diff) > (se_diff * 2))
    } else (
      tibble(compare = NA, 
             models = NA)  
    )
  }) %>% 
  ungroup() %>% 
  select(test_subset, symp, grouping, model_version, distribution, elpd_diff, se_diff, rank, significant, weight) 

saveRDS(loo_compare, file = here("generated_data/all_loo_compared.rds"))

library(kableExtra)
loo_compare %>% 
  kable(align = "c") %>%
  kable_styling(full_width = F) %>%
  column_spec(1, bold = T) %>%
  collapse_rows(columns = 1:2, valign = "top")


# Stacked summaries -------------------------------------------------------
# Compute model stack summaries for quantities of interest

# Trajectories of CDF
cdf_stack <- model_stack(post_draws = inc_period_q_traj,
                         loo_compare = loo_compare,
                         var_cols = c("ps", "grouping", "group"),
                         n_draws = 300)

saveRDS(cdf_stack, file = "generated_data/quantile_traj_stack.rds")


# Quantiles of interest
quantile_stack_summary <- model_stack_summary(post_draws = inc_period_q_draws,
                                              loo_compare = loo_compare,
                                              var_cols = c("ps", "grouping", "group"),
                                              n_draws = 4000)

saveRDS(quantile_stack_summary, file = "generated_data/quantile_stack_summary.rds")

pd <- position_dodge(.05)

quantile_stack_summary %>% 
  filter(ps <= 0.95, ps >= 0.05) %>% 
  ggplot(aes(x = ps, y = mean, color = symp)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), position = pd) +
  facet_wrap(~test_subset) +
  theme_bw() +
  coord_flip()

# Probability of exceedence
p_exceed_stack_summary <- model_stack_summary(post_draws = prop_exceed_draws,
                                              loo_compare = loo_compare,
                                              var_cols = c("thresh", "grouping", "group"),
                                              n_draws = 4000)

saveRDS(p_exceed_stack_summary, file = "generated_data/p_exceed_stack_summary.rds")

