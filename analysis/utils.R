# This script contains functions used in the analysis


# Data pulling ------------------------------------------------------------

#' read_raw_data
#' Read raw csv files stored in onedrive folder
#'
#' @param data_paths paths to onedrive folder as returned by `set_paths`
#' @param file_name file to read
#'
#' @return tibble with data
#'
read_raw_data <- function(raw_data_path = here::here("data/raw_data"),
                          file_name = NULL) {
  
  if (is.null(raw_data_path)) {
    # If missing set path to onedrive data folder
    raw_data_paths <- set_paths()$raw_data
  }
  
  # Read in raw file
  dat <- readr::read_csv(here::here(raw_data_path, file_name))
  
  dat
}


#' simplify_colnames
#' remove questionnaire groups from column names
#' 
#' @param df 
#'
#' @return modified dataframe with simpler column names
#' 
simplify_colnames <- function(df) {
  df %>% 
    magrittr::set_colnames(str_remove(colnames(df), "_f[0-9]+"))
}

# Data processing ---------------------------------------------------------

#' get_age_cuts
#' 
get_age_cuts <- function() {
  c(0, 1, 5, 15, Inf)
}

#' get_age_cat_levels
#'
get_age_cat_levels <- function() {
  c("[0,1)" = "Under 1 year", "[1,5)" = "1-4 years", 
    "[5,15)" = "5-14 years", "[15,Inf]" = "≥ 15 years")
}

get_sex_levels <- function() {
  c("female", "male")
}

#' add_age_num
#' Adds a column with numeric age in years
#' 
#' @param df dataframe to modify
#' @param col_age name of the numeric column containing age
#' @param col_age_unit vector of age cuts to apply
#'
#' @return
#'
add_age_num <- function(df,
                        col_age = "age",
                        col_age_unit = "age_unit") {
  
  df %>% 
    dplyr::mutate(
      age_num = case_when(
        get(col_age_unit) == "month" ~ get(col_age)/12,
        T ~ !!rlang::sym(col_age))
    )
}

#' add_age_cat
#' Sets age categories in a dataframe into a new column called `age_cat`
#' 
#' @param df dataframe to modify
#' @param col_age name of the numeric column containing age
#' @param age_cuts vector of age cuts to apply
#'
#' @return
#'
add_age_cat <- function(df,
                        col_age = "age_num",
                        age_cuts = get_age_cuts()) {
  
  df %>% 
    dplyr::mutate(age_cat = cut(!!rlang::sym(col_age), age_cuts, 
                                include.lowest = T, right = F) %>% 
                    factor(levels = names(get_age_cat_levels()),
                           labels = get_age_cat_levels()))
  
}


#' extract_numeric
#' extract numeric values from string string variable
#' 
#' @param string variable to extract values from
#'
#' @return numeric values
#'

extract_numeric <- function(string = "gps") {
  
  # extract numeric values 
  
  values <- str_extract_all(string, "-?\\d+\\.\\d+")[[1]]
  values <- as.numeric(values)
  
  return(values)
}


#' get_gps
#' create longitude and latitude variables from values contained in a string variable
#' 
#' @param df data frame to modify
#' @param col_gps name of the character column containing gps coordinates
#'
#' @return data frame with longitude and latitude columns
#'

get_gps <- function(df, col_gps = "gps") {
  
  df <- df %>%
    rowwise() %>%
    mutate(latitude = extract_numeric(col_gps)[1],
           longitude = extract_numeric(col_gps)[2]) %>%
    ungroup() 
  
  return(df)
}



#' colname_parser_checks
#' Check functions for parsing colnames (for use in pivot_longer_parts)
#' 
#' @param col_names colanems to parse
#' @param parsed_names parsed colnames (as returned by parser)
#'
colname_parser_checks <- function(col_names, parsed_names) {
  collapsed <- str_c(parsed_names, collapse = ", ")
  
  if (any(duplicated(parsed_names))) {
    stop("Colname parsing caused duplicates, please check: ", collapsed)
  }
  
  if (length(parsed_names) != length(col_names)) {
    stop("Colname parsing dimension missmatch, please check: ", collapsed)
  }
}

#' genexp_test_colname_parser
#' Colname parser for genexpert test fields
#' 
#' @param col_names 
#'
genexp_test_colname_parser <- function(col_names) {
  parsed_names <- str_c("test_", str_extract(col_names, "result$|^dt"))
  parsed_names[is.na(parsed_names)] <- "test_bool"
  
  # Run checks
  colname_parser_checks(col_names, parsed_names)
  
  parsed_names
}

#' pivot_longer_parts
#' Pivot dataframe longer by grouping columns by pattern
#' 
#' @param df dataframe to pivot
#' @param keys columns to use as keys while pivoting
#' @param patterns unique string patterns to use for grouping columns
#' @param colname_parser colname parser
#'
pivot_longer_parts <- function(df, 
                               keys,
                               patterns,
                               colname_parser) {
  
  map_dfr(patterns, function(x) {
    
    # Get keys
    keys_df <- df %>% select(all_of(keys))
    
    # Get data to pivot
    dat_df <- df %>% 
      select(contains(x)) %>% 
      magrittr::set_colnames(colname_parser(colnames(.)))
    
    # Bind and return
    bind_cols(keys_df, dat_df) %>% 
      mutate(what = x)
  })
}


#' add_epi_week
#' Add epidemiologic week
#' 
#' @param df 
#' @param date_col 
#'
#' @return
#' @export
#'
#' @examples
add_epi_week <- function(df, date_col = "date") {
  
  df %>% 
    mutate(week = lubridate::epiweek(!!rlang::sym(date_col)),
           year = lubridate::epiyear(!!rlang::sym(date_col)),
           epiweek = str_c(year, week, sep = "-"),
           epiweek_date = as.Date(str_c(epiweek, "-1"), format = "%Y-%U-%u"))
}

get_test_options <- function() {
  c("all", "positive", "negative")
}

get_symptom_options <- function() {
  c("rash", "fever", "any")
}

#' make_dic_data
#' format data for fit with coarseDataTools
#' @param df 
#' @param from_col 
#' @param to_col 
#'
#' @return
#' @export
#'
#' @examples
make_dic_data <- function(df, 
                          from_col = "from_date", 
                          to_col = "to_date",
                          min_exp_date) {
  
  df %>% 
    mutate(from_date = case_when(is.na(!!rlang::sym(from_col)) ~ min_exp_date,
                                 !!rlang::sym(to_col) < !!rlang::sym(from_col) ~ min_exp_date,
                                 T ~ !!rlang::sym(from_col))) %>% 
    mutate(SL = !!rlang::sym(to_col),
           # no uncertainty in symptom onset
           SR = !!rlang::sym(to_col),  
           # right censored observation
           EL = case_when(contact_occurence == "once" ~ from_date,
                          T ~ min_exp_date),
           ER = case_when(is.na(contact_occurence) ~ !!rlang::sym(to_col),
                          from_date == min_exp_date ~ !!rlang::sym(to_col),
                          T ~ from_date),  
           # Transform to numeric values
           across(c("SL", "SR", "EL", "ER"), ~ as.numeric(difftime(., min_exp_date, units = "days"))),
           # Add half a day when symptom onset reported on same day as exposure time
           epsilon = case_when(SL == ER ~ 0.5,
                               T ~ 0),
           SL = SL + epsilon,
           SR = SR + epsilon,
           E_int = ER-EL,
           S_int = SR-SL,
           max_inc_int = SR-EL,
           min_inc_int = SL-ER) 
}

diff_time_days <- function(x, y) {
  difftime(x, y, units = "days") %>% as.numeric()
}

#' tidy_register_id
#' standardize register id as a 3-character string with 0 padding
#' @param df 
#'
tidy_register_id <- function(df) {
  df %>% 
    mutate(register_id = str_pad(as.character(as.numeric(register_id)), 
                                 side = "left", pad = "0", width = 4))
}

#' tidy_contact_labels
#' tidy duplicated labels with spelling variation for suspected contacts 
#' 
#' @param dat clinical contact dataframe
#'
tidy_contact_labels <- function(dat){
  
  # voisin
  dat$other_relationship[dat$other_relationship %in% 
                           c("Voisin","voisine","Voisine")] <- "voisin"
  
  # enfants du voisin
  dat$other_relationship[dat$other_relationship %in% 
                           c("enfant du voisin","3 enfants du voisin")] <- "enfants du voisin"
  
  # "Non prC)cisC)" 
  dat$other_relationship[dat$other_relationship=="Non prC)cisC)"] <- "Non precise"
  
  return(dat)
}

# External data -----------------------------------------------------------

get_lesion_resolution_pars <- function() {
  # Back-infer log-normal time to lesion resolution
  # Values from https://www.medrxiv.org/content/10.1101/2024.11.18.24316975v1.full.pdf
  objfun <- function(log_pars, quantiles, values) {
    sim_values <- qlnorm(1-quantiles, meanlog = log_pars[1], sdlog = exp(log_pars[2]))
    sum((values - sim_values)^2)
  }
  
  set.seed(32423)
  best_pars <- optim(par = c(log(10), log(1)),
                     fn = objfun,
                     quantiles = c(0.95, 0.75, .5, .25),
                     values = c(5, 11, 14, 20))
  
  best_pars
}

# Checks ------------------------------------------------------------------

#' run_checks_clinical_data
#' checks for clinical data
#' @param dat 
#'
run_checks_clinical_data <- function(dat) {
  
  # Age categories
  expect_false(any(is.na(dat$age_cat)))
  
  # Unique register ids (being fixed in https://github.com/HopkinsIDD/mpox-uvira-sprint/issues/3#issue-2699504737)
  # expect_false(any(duplicated(dat$id)))
  
  # Dates
  expect_gte(min(dat$dt_admission, na.rm = T), as.Date("2024-05-01"))
  expect_lte(max(dat$dt_admission, na.rm = T), as.Date("2024-12-31"))
  
  cat("---- Passed clinical data checks. \n")
}

#' run_checks_test_data
#' checks for genexpert test data
#' @param dat 
#'
run_checks_test_data <- function(dat) {
  cat("---- No test data checks yet. \n")
  
}

#' run_checks_contact_data
#' checks for contact data
#' @param dat 
#'
run_checks_contact_data <- function(dat) {
  cat("---- No contact data checks yet. \n")
  
}


#' run_checks_full_data
#' checks for combiend clinical, test, and contact data
#' @param dat 
#'
run_checks_full_data <- function(dat) {
  cat("---- No full data checks yet. \n")
  
}

#' run_checks_options
#'
#' @param opt 
#'
#' @return
#' @export
#'
#' @examples
run_checks_options <- function(opt) {
  
  # Check options
  if (!(opt$test_subset %in% get_test_options())) {
    stop("Please set test subset in: ", str_c(get_test_options(), collapse = ", "))
  }
  
  if (!(opt$symptom_type %in% get_symptom_options())) {
    stop("Please set symptom type to: ", str_c(get_symptom_options(), collapse = ", "))
  }
  
}

# Prep for stan utils -----------------------------------------------------

make_starts <- function(map_to_i) {
  u_i <- unique(map_to_i)
  starts <- rep(0, length(u_i))
  for (i in 1:length(u_i)) {
    starts[i] <- which(map_to_i == i)[1]
  }
  starts
}

make_ends <- function(map_to_i) {
  u_i <- unique(map_to_i)
  ends <- rep(0, length(u_i))
  for (i in 1:length(u_i)) {
    ends[i] <- which(map_to_i == i) %>% last()
  }
  ends
}

make_map <- function(vec, u_ref) {
  map_dbl(vec, ~ which(u_ref == .))
}


# File names --------------------------------------------------------------

#' make_filename
#'
#' @param basename 
#' @param outdir 
#' @param opt 
#' @param parse_opt 
#' @param type 
#' @param verbose 
#'
make_filename <- function(basename, 
                          outdir = here("generated_data/"), 
                          opt = NULL, 
                          parse_opt = NULL,
                          type = "rds",
                          verbose = TRUE) {
  
  if (!is.null(parse_opt)) {
    parsed_opt <- str_c("_", parse_opt(opt))
  } else {
    parsed_opt <- ""
  }
  
  file_name <- str_glue("{outdir}{basename}{parsed_opt}.{type}")
  
  if (verbose) {
    cat("-- file name:", file_name, "\n")
  }
  
  file_name
}

#' make_clinical_filename
#' filename for processed clinical data
#' @param outdir 
#'
make_clinical_filename <- function(outdir = here("generated_data/")) {
  make_filename(basename = "processed_clinical_data",
                outdir = outdir,
                type = "rds")
}

#' make_contacts_filename
#' filename for processed contact data
#' @param outdir 
#'
make_contacts_filename <- function(outdir = here("generated_data/")) {
  make_filename(basename = "processed_contacts_data",
                outdir = outdir,
                type = "rds")
}

#' make_test_data_filename
#' filename for processed test data
#' @param outdir 
#'
make_test_data_filename <- function(outdir = here("generated_data/")) {
  make_filename(basename = "processed_test_data",
                outdir = outdir,
                type = "rds")
}

#' make_test_filename
#' filename for combined test, clinical and contact data
#' @param outdir 
#'
make_combined_filename <- function(outdir = here("generated_data/")) {
  make_filename(basename = "processed_combined_data",
                outdir = outdir,
                type = "rds")
}


#' parse_stan_data_opt
#' parser of option for options related to running the stan model for incubation period
#' @param opt 
#'
parse_inc_period_stan_opt <- function(opt) {
  res <- str_glue("usedata-{opt$use_data}_symp-{opt$symptom_type}_test-{opt$test_subset}_group-{opt$group_column}_mversion-{opt$model_version}_dist-{opt$dist_type}")
  
  if (is.null(res)) {
    stop("Failed parsing options.")
  }
  
  res
}

#' parse_stan_data_opt
#' parser of option for options related to running the stan model for incubation period
#' @param opt 
#'
parse_inc_period_dic_stan_opt <- function(opt) {
  res <- str_glue("usedata-{opt$use_data}_symp-{opt$symptom_type}_test-{opt$test_subset}_group-{opt$group_column}_mversion-{opt$model_version}_dist-{opt$dist_type}")
  
  if (is.null(res)) {
    stop("Failed parsing options.")
  }
  
  res
}

#' make_inc_period_stan_data_filename
#' filename for stan data for incubation period
#' 
#' @param outdir 
#' @param opt 
#'
make_inc_period_stan_data_filename <- function(outdir = here("generated_data/"),
                                               opt) {
  make_filename(basename = "inc_period_stan_data",
                outdir = outdir,
                opt = opt,
                parse_opt = parse_inc_period_stan_opt,
                type = "rds")
}

#' make_inc_period_dic_data_filename
#' filename for stan data for incubation period
#' 
#' @param outdir 
#' @param opt 
#'
make_inc_period_dic_data_filename <- function(outdir = here("generated_data/"),
                                              opt) {
  make_filename(basename = "inc_period_dic_data",
                outdir = outdir,
                opt = opt,
                parse_opt = parse_inc_period_stan_opt,
                type = "rds")
}

#' make_inc_period_stan_data_filename
#' filename for stan fit for incubation period
#' 
#' @param outdir 
#' @param opt 
#'
make_inc_period_stan_fit_filename <- function(outdir = here("generated_data/"),
                                              opt) {
  
  make_filename(basename = "inc_period_stan_fit",
                outdir = outdir,
                opt = opt,
                parse_opt = parse_inc_period_stan_opt,
                type = "rds")
}

#' make_inc_period_dic_fit_filename
#' filename for stan fit for incubation period
#' 
#' @param outdir 
#' @param opt 
#'
make_inc_period_dic_fit_filename <- function(outdir = here("generated_data/"),
                                             opt) {
  
  make_filename(basename = "inc_period_dic_fit",
                outdir = outdir,
                opt = opt,
                parse_opt = parse_inc_period_stan_opt,
                type = "rds")
}

#' make_inc_period_genquant_filename
#' filename for stan generated quantities for incubation period
#' 
#' @param outdir 
#' @param opt 
#'
make_inc_period_genquant_filename <- function(outdir = here("generated_data/"),
                                              opt) {
  
  make_filename(basename = "inc_period_genquant",
                outdir = outdir,
                opt = opt,
                parse_opt = parse_inc_period_stan_opt,
                type = "rds")
}


parse_inc_period_stan_filename <- function(str) {
  tibble(
    use_data = str_extract(str, "(?<=usedata-)[TRUE|FLASE]") %>% as.logical(),
    test_subset = str_extract(str, "(?<=test-)[a-z]+"),
    # sd_prior = str_extract(str, "(?<=sdprior-)(.)*(?=_)") %>% as.numeric(),
    symp = str_extract(str, "(?<=symp-)[a-z]+"),
    grouping = str_extract(str, "(?<=group-)(.)*(?=_mversion)"),
    model_version = str_extract(str, "(?<=mversion-)(.)*(?=_dist)"),
    distribution = str_extract(str, "(?<=dist-)(.)*(?=\\.)")
  )
}

# Plotting ----------------------------------------------------------------


plot_posteriors <- function(fit_files) {
  draws <- map_df(fit_files, function(x) {
    res <- readRDS(x)
    stan_data <- readRDS(str_replace(x, "fit", "data"))
    
    res$draws(c("mu", "sigma")) %>%
      as_draws_df() %>% 
      as_tibble() %>% 
      pivot_longer(contains("["),
                   names_to = "variable",
                   values_to = "value") %>% 
      mutate(var = str_extract(variable, "mu|sigma"),
             group = stan_data$u_groups[as.numeric(str_extract(variable, "[0-9]+"))]) %>% 
      select(-variable) %>% 
      bind_cols(parse_inc_period_stan_filename(x))
  })
  
  
  draws %>% 
    filter(use_data, var == "mu") %>% 
    ggplot(aes(x = value)) + 
    geom_density(aes(color = test_subset)) +
    geom_density(data = draws %>% 
                   filter(!use_data, var == "mu") %>% 
                   mutate(test_subset = ifelse(!use_data, "prior", test_subset)) %>% 
                   select(-symp)) +
    facet_grid(sd_prior ~ symp, scales = "free", labeller = label_both) +
    theme_bw() +
    scale_color_manual(values = c("blue", "red", "black")) +
    coord_cartesian(xlim = c(0, 20)) +
    labs(x = "mean incubation period")
}


table_posteriors <- function(fit_files) {
  
  summaries <- map_df(fit_files, function(x) {
    res <- readRDS(x)
    stan_data <- readRDS(str_replace(x, "fit", "data"))
    
    res$summary(c("mu", "sigma"), custom_summaries()) %>%
      mutate(var = str_extract(variable, "mu|sigma"),
             group = stan_data$u_groups[as.numeric(str_extract(variable, "[0-9]+"))]) %>% 
      select(-variable) %>% 
      bind_cols(parse_inc_period_stan_filename(x))
  }) %>% 
    rename(lo = q2.5, hi = q97.5) %>% 
    make_estimate_text() %>% 
    mutate(test_subset = ifelse(!use_data, "prior", test_subset)) %>% 
    select(var, test_subset, sd_prior, txt) %>% 
    pivot_wider(names_from = "var",
                values_from = "txt") %>% 
    arrange(test_subset, sd_prior) 
  
  summaries %>% 
    select(sd_prior, mu, sigma) %>% 
    kable(col.names = c("sd of prior", "mean (95 CrI)", "sd (95 CrI)")) %>% 
    kable_styling(full_width = F) %>% 
    pack_rows(index = table(summaries$test_subset))
}


# Post processing --------------------------------------------------------

format_number <- function(x) {
  formatC(x, format = "f", digits = 1, big.interval = 3, big.mark = ",")
}

format_int <- function(x) {
  formatC(x, format = "f", digits = 0, big.interval = 3, big.mark = ",")
}

make_estimate_text <- function(df,
                               pct = FALSE) {
  
  # multiplier if the estimate is a percent
  multi <- ifelse(pct, 100, 1)
  
  # symbol
  s <- ifelse(pct, "% ", " ")
  
  df %>% 
    mutate(txt = str_glue("{format_number(mean*multi)}{s}({format_number(lo*multi)}-{format_number(hi*multi)})"))
}


#' custom_summaries
#' Custom summaries to get the 95% CrI
#'
custom_summaries <- function() {
  
  c(
    "mean", "median", "custom_quantile2",
    posterior::default_convergence_measures(),
    posterior::default_mcse_measures()
  )
}

#' cri_interval
#' The Credible interval to report in summaries
#' 
cri_interval <- function() {
  c(0.025, 0.975)
}

#' custom_quantile2
#' quantile functoin with custom cri
#'
#' @param x 
#' @param cri 
#'
#'
custom_quantile2 <<- function(x, cri = cri_interval()) {
  posterior::quantile2(x, probs = cri)
}


#' get_dic_params
#' get location and scale from coarseDataTools dic.fit object
#' @param dic_fit 
#'
#' @return
#' @export
#'
#' @examples
get_dic_params <- function(dic_fit) {
  dic_fit@ests[c(1:2), -4] %>% 
    as_tibble() %>% 
    rename(mean = est, q2.5 = CIlow, q97.5 = CIhigh)
}

#' get_dic_quantiles
#' get quantiles from coarseDataTools dic.fit object
#' @param dic_fit 
#'
#' @return
#' @export
#'
#' @examples
get_dic_quantiles <- function(dic_fit) {
  dic_fit@ests[-c(1:2), -4] %>% 
    as_tibble() %>% 
    mutate(p = rownames(dic_fit@ests)[-c(1:2)] %>% 
             str_remove("p") %>% 
             as.numeric() %>% 
             {./100}) %>% 
    rename(mean = est, q2.5 = CIlow, q97.5 = CIhigh) 
}


#' compute_quantiles
#' Compute the quantiles for a target distribution
#' 
#' @param stan_fit 
#' @param ps 
#' @param dist 
#' @param u_groups
#'
compute_quantiles <- function(stan_fit, 
                              ps, 
                              dist, 
                              u_groups,
                              subset = FALSE) {
  
  if (dist == 0) {
    f <- qlnorm
  } else if (dist == 1){
    f <- qgamma
  } else if (dist == 2) {
    f <- qweibull
  }
  
  stan_fit$draws("pars") %>% 
    as_draws() %>% 
    as_draws_df() %>% 
    as_tibble() %>% 
    pivot_longer(cols = contains("pars"),
                 names_to = "variable") %>% 
    mutate(what = str_c("par", str_extract(variable, "(?<=\\[)[0-9]+")),
           g = str_extract(variable, "(?<=\\,)[0-9]+") %>% as.numeric(),
           group = u_groups[g]) %>% 
    {
      dat <- .
      if (subset) {
        filter(dat, .draw %in% sample(1:4000, 100))
      } else {
        dat
      }
    } %>% 
    select(-variable, -.chain, -.iteration) %>% 
    pivot_wider(names_from = "what",
                values_from = "value") %>% 
    group_by(group, .draw) %>% 
    group_modify(
      function(x, y){
        tibble(
          ps = ps,
          value = f(ps, x$par1, x$par2)
        )
      }
    )
}

#' summarise_quantiles
#'
#' @param stan_fit 
#' @param ps 
#' @param dist 
#' @param u_groups 
#'
summarise_quantiles <- function(stan_fit, 
                                ps, 
                                dist, 
                                u_groups) {
  
  compute_quantiles(stan_fit = stan_fit, 
                    ps = ps, 
                    dist = dist, 
                    u_groups = u_groups) %>% 
    group_by(group, ps) %>% 
    summarise_draws() %>% 
    ungroup()
}

#' summarise_draws
#' helper to summarise draws when in tibble format
#' @param df 
#'
summarise_draws <- function(df) {
  df %>% 
    summarise(
      mean = mean(value, na.rm = T),
      median = median(value, na.rm = T),
      q2.5 = quantile(value, 0.025, na.rm = T),
      q97.5 = quantile(value, 0.975, na.rm = T),
      .groups = "keep")
}

#' model_stack
#' Compute model stack summaries
#' 
#' @param post_draws 
#' @param loo_compare 
#' @param n_draws 
#'
model_stack <- function(post_draws, 
                        loo_compare,
                        n_draws = 4000,
                        var_cols = NULL) {
  
  group_cols <- c("test_subset", "symp", var_cols)
  
  inner_join(post_draws, loo_compare) %>% 
    group_by_at(group_cols) %>% 
    group_modify(function(x, y) {
      x %>% 
        sample_n(size = n_draws, 
                 weight = weight, 
                 replace = F)
    })
}

#' model_stack_summary
#' Compute model stack summaries
#' 
#' @param post_draws 
#' @param loo_compare 
#' @param n_draws 
#'
model_stack_summary <- function(post_draws, 
                                loo_compare,
                                n_draws = 4000,
                                var_cols = NULL) {
  
  group_cols <- c("test_subset", "symp", var_cols)
  
  model_stack(post_draws = post_draws, 
              loo_compare = loo_compare,
              n_draws = n_draws,
              var_cols = var_cols) %>% 
    group_by_at(group_cols) %>% 
    summarise_draws() %>% 
    ungroup()
}

#' add_factors_to_results
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
add_factors_to_results <- function(df) {
  df %>% 
    mutate(
      grouping = factor(grouping, 
                        levels = c("none", "age_cat_simple", "sex", "period", "hospitalized", "exposure_cat_simple", "exposure_relation_simple"),
                        labels = c("overall", "age group", "sex", "time\nperiod", "hospitalization", "reported\nexposure type", "reported\nrelation type")),
      group = factor(group,
                     levels = c("all", "u15", "15+", "female", "male", "early", "late", "yes", "no", "sexual", "nonsexual", "sexual_rel", "nonsexual_rel"),
                     labels = c("overall", "under 15", "15+", "female", "male",
                                "early\n(before Sep 1st)", "late\n(after Sep 1st)",
                                "inpatient", "outpatient", "sexual", "non-sexual", "sexual relation", "non-sexual relation")),
      test_subset = factor(test_subset, 
                           levels = c("positive", "all"),
                           labels = c("tested positive", "all patients")),
      distribution_label = factor(distribution, 
                                  levels = c("L", "G", "W"),
                                  labels = c("log-normal", "gamma", "weibull"))
    )
}

#' map_df_progress
#' https://www.jamesatkins.com/posts/progress-bar-in-purrr-map-df/
#' @param .x 
#' @param .f 
#' @param ... 
#' @param .id 
#'
map_df_progress <- function(.x, .f, ..., .id = NULL) {
  .f <- purrr::as_mapper(.f, ...)
  pb <- progress::progress_bar$new(total = length(.x), force = TRUE)
  
  f <- function(...) {
    pb$tick()
    .f(...)
  }
  purrr::map_df(.x, f, ..., .id = .id)
}



#' map_postprocess
#' generic function to apply postprocessing function to list of stan fit files
#' 
#' @param files 
#' @param fun 
#'
map_postprocess <- function(files, fun, .options = furrr::furrr_options(globals = TRUE)) {
  future_map_dfr(files, function(x) {
    stan_fit <- readRDS(x)
    stan_data <- readRDS(str_replace(x, "stan_fit|genquant", "stan_data"))
    
    fun(stan_fit, stan_data) %>% 
      bind_cols(parse_inc_period_stan_filename(x))
  },
  .options = .options)
}

#' postprocess_mu_sigma
#' extract mean and sd of distribution
#' @param stan_fit 
#' @param stan_data 
#' 
postprocess_mu_sigma <- function(stan_fit, stan_data) {
  stan_fit$summary(c("mu", "sigma"), custom_summaries()) %>%
    mutate(var = str_extract(variable, "mu|sigma"),
           group = stan_data$u_groups[as.numeric(str_extract(variable, "[0-9]+"))]) %>% 
    select(-variable)
}

#' postprocess_quantiles_traj
#' compute target quantiles
#' @param stan_fit 
#' @param stan_data 
#'
postprocess_quantiles_traj <- function(stan_fit, stan_data) {
  compute_quantiles(stan_fit = stan_fit,
                    ps = seq(0.01, 0.99, length.out = 50),
                    dist = stan_data$stan_data$dist,
                    u_groups = stan_data$u_groups,
                    subset = TRUE)
}

#' postprocess_quantiles_target_draws
#' compute quantiles for model stacking
#' @param stan_fit 
#' @param stan_data 
#'
postprocess_quantiles_target_draws <- function(stan_fit, 
                                               stan_data) {
  compute_quantiles(stan_fit = stan_fit,
                    ps = stan_data$stan_data$p_target,
                    dist = stan_data$stan_data$dist,
                    u_groups = stan_data$u_groups,
                    subset = FALSE)
}



#' postprocess_quantiles
#' compute target quantiles
#' @param stan_fit 
#' @param stan_data 
#'
postprocess_quantiles <- function(stan_fit, stan_data) {
  summarise_quantiles(stan_fit = stan_fit,
                      ps = stan_data$stan_data$p_target,
                      dist = stan_data$stan_data$dist,
                      u_groups = stan_data$u_groups) 
}


#' postprocess_cdfs
#' compute incubation period cdfs
#' 
#' @param stan_fit 
#' @param stan_data 
#'
postprocess_cdfs <- function(stan_fit, stan_data) {
  summarise_quantiles(stan_fit = stan_fit,
                      ps =  seq(0.01, 0.99, length.out = 50),
                      dist = stan_data$stan_data$dist,
                      u_groups = stan_data$u_groups) 
}

#' postprocess_loo
#'
#' @param stan_fit 
#' @param fun 
#'
postprocess_loo <- function(stan_fit, fun) {
  tibble(loo = list(loo::loo(stan_fit$draws("log_lik") %>% 
                               as_draws() %>% 
                               subset_draws(draw = 1:3000)
                             )))
}


#' postprocess_prop_exceed
#' extract mean and sd of distribution
#' @param stan_fit 
#' @param stan_data 
#' 
postprocess_prop_exceed <- function(stan_fit, stan_data) {
  threshs <- factor(c(14, 21, 28))
  stan_fit$summary("prop_exceed", custom_summaries()) %>%
    mutate(thresh = threshs[as.numeric(str_extract(variable, "(?<=\\,)[0-9]+"))],
           group = stan_data$u_groups[as.numeric(str_extract(variable, "(?<=\\[)[0-9]+"))]) %>% 
    select(-variable)
}

#' postprocess_prop_exceed_draws
#' extract mean and sd of distribution
#' @param stan_fit 
#' @param stan_data 
#' 
postprocess_prop_exceed_draws <- function(stan_fit, stan_data) {
  threshs <- factor(c(14, 21, 28))
  stan_fit$draws("prop_exceed") %>%
    posterior::as_draws() %>% 
    posterior::as_draws_df() %>% 
    as_tibble() %>% 
    pivot_longer(cols = contains("p_"),
                 names_to = "variable") %>% 
    mutate(thresh = threshs[as.numeric(str_extract(variable, "(?<=\\,)[0-9]+"))],
           group = stan_data$u_groups[as.numeric(str_extract(variable, "(?<=\\[)[0-9]+"))]) %>% 
    select(-variable)
}


#' kable_res
#' Function to make a simple tidy kable
#'
#'
#' @param df 
#' @param var_col 
#' @param lo_col 
#' @param hi_col 
#' @param pct 
#'
kable_res <- function(df, 
                      var_col, 
                      lo_col = "q2.5", 
                      hi_col = "q97.5", 
                      pct = FALSE) {
  df2 <- df %>% 
    rename(lo = !!rlang::sym(lo_col), hi = !!rlang::sym(hi_col)) %>% 
    make_estimate_text(pct = pct) %>% 
    select(test_subset, symp, model_version, grouping, group, distribution_label, txt, !!rlang::sym(var_col)) %>% 
    pivot_wider(names_from = var_col,
                values_from = "txt") %>% 
    arrange(test_subset, symp, model_version, grouping, group, distribution_label)
  
  df2 %>% 
    kable(align = "c") %>%
    kable_styling(full_width = F) %>%
    column_spec(1, bold = T) %>%
    collapse_rows(columns = 1:5, valign = "top")
}


# Stan helpers ------------------------------------------------------------


make_stan_init <- function(stan_data, v){
  
  J <- stan_data$J
  J_incid <- ifelse(J > 1 && stan_data$incid_by_group == 1, J, 1)
  multigroup <- J > 1
  P <- stan_data$P
  d <- stan_data$dist + 1
  
  bounds <- list(
    sigma = list(
      left = c("L" = 0.5, "G" = 5, "W" = 4),
      right = c("L" = 0.8, "G" = 10, "W" = 7)
    )
  )
  
  if (v == "v1") {
    log_rho <- matrix(rnorm(P * J_incid, 0, .1), nrow = P, ncol = J_incid)
    
    for (i in 1:J_incid) {
      log_rho[, i] <- rnorm(P, seq(-2, 2, length.out = P), .1)
    }
    
    init <- map(1:4, ~ list(
      log_mu_std = rnorm(J, 0, .3),
      sigma = runif(J, bounds$sigma$left[d], bounds$sigma$right[d]),
      log_mu_mu = rnorm(multigroup, log(9), .1),
      sd_mu = runif(multigroup, 0.5, 1),
      mu_logit_lambda = rnorm(J, 2, .3),
      log_rho = log_rho,
      mu_rho = rnorm(J_incid, -1, .3),
      sd_rho = runif(J_incid, .5, 1),
      log_mu_reporting = rnorm(J, log(5), .1),
      sigma_reporting = runif(J, .7, .9),
      inv_od = runif(1, 0.1, 0.5),
      logit_phi = runif(1, -1, -0.5)
    ))
  } else {
    init <- map(1:4, ~ list(
      log_mu_std = rnorm(J, 0, .1),
      sigma = runif(J, bounds$sigma$left[d], bounds$sigma$right[d]),
      log_mu_mu = rnorm(multigroup, log(9), .1),
      sd_mu = runif(multigroup, 0.5, 1),
      logit_phi = runif(1, -1, -0.5)
    ))
  }
  
  init
}

