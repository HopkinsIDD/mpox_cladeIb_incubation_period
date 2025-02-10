# This script runs the stan model


# Preamble ---------------------------------------------------------------
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(purrr)
library(stringr)
library(optparse)
library(here)
library(cmdstanr)
library(posterior)
library(bayesplot)

# Functions for analysis
source("analysis/utils.R")

# User-supplied options (if any)
option_list <- list(
  make_option(c("-c", "--config"), 
              default = NULL, action ="store", type = "character", 
              help = "Do plots while processing data."),
  make_option(c("-y", "--symptom_type"), 
              default = "rash", action ="store", type = "character", 
              help = "type of symptom to consider"),
  make_option(c("-x", "--group_column"), 
              default = "none", action ="store", type = "character", 
              help = "Columns to use for grouping"),
  make_option(c("-t", "--test_subset"), 
              default = "positive", action ="store", type = "character", 
              help = "Subset of genexpert test resutls to use"),
  make_option(c("-v", "--model_version"), 
              default = "v0", action ="store", type = "character", 
              help = "Stan model version"),
  make_option(c("-d", "--use_data"), 
              default = TRUE, action ="store", type = "logical", 
              help = "Whether to use data in inference or only prior"),
  make_option(c("-r", "--run_stan"), 
              default = TRUE, action ="store", type = "logical", 
              help = "F to re-run stan model"),
  make_option(c("-g", "--run_genquant"), 
              default = TRUE, action ="store", type = "logical", 
              help = "Whether to make generated quantities"),
  make_option(c("-z", "--dist_type"), 
              default = "W", action ="store", type = "character", 
              help = "Distribution type from incubation period. One of L, G, W")
)

# Parse options
opt_ <- parse_args(OptionParser(option_list = option_list))

if (!is.null(opt_$config)) {
  opt <- yaml::read_yaml(opt_$config)
} else {
  opt <- opt_
}

opt$run_genquant <- opt_$run_genquant
opt$run_stan <- opt_$run_stan

# Run checks and print for logs
run_checks_options(opt)
print(opt)

# Load data ---------------------------------------------------------------
stan_data <- readRDS(here(make_inc_period_stan_data_filename(opt = opt)))

# Run stan model ----------------------------------------------------------
if (opt$run_stan | !file.exists(make_inc_period_stan_fit_filename(opt = opt))) {
  
  stan_model <- cmdstan_model(here(str_glue("analysis/stan/infer_incubation_period_{opt$model_version}_trunc2.stan")))
  
  # Make starting points for sampler
  init <- make_stan_init(stan_data = stan_data$stan_data,
                         v = opt$model_version)
  
  # Draw samples
  stan_fit <- stan_model$sample(data = stan_data$stan_data,
                                init = init,
                                chains = 4,
                                parallel_chains = 4,
                                iter_warmup = 250,
                                iter_sampling = 1000,
                                refresh = 50, 
                                save_warmup = FALSE)
  
  # Save stan fit
  stan_fit$save_object(make_inc_period_stan_fit_filename(opt = opt))
  
} else {
  # Load
  stan_fit <- readRDS(make_inc_period_stan_fit_filename(opt = opt))
}

# Run genquant ------------------------------------------------------------
if (opt$run_genquant) {
  
  gen_model <- cmdstan_model(here(str_glue("analysis/stan/infer_incubation_period_{opt$model_version}_generate.stan")))
  
  # Run generated quantities
  genquant <- gen_model$generate_quantities(fitted_params = stan_fit,
                                            data = stan_data$stan_data, 
                                            parallel_chains = 4)
  
  # Save generated quantities
  genquant$save_object(make_inc_period_genquant_filename(opt = opt))
  
}
