# This script prepares the data for incubation period estimation


# Preamble ---------------------------------------------------------------

# Load packages
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)
library(forcats)
library(testthat)
library(optparse)
library(here)
# library(coarseDataTools)

# Functions for analysis
source("analysis/utils.R")

# User-supplied options (if any)
option_list <- list(
  make_option(c("-c", "--config"), 
              default = NULL, action ="store", type = "character", 
              help = "Do plots while processing data."),
  make_option(c("-p", "--do_plots"), 
              default = T, action ="store", type = "logical", 
              help = "Do plots while processing data."),
  make_option(c("-r", "--redo_data_processing"), 
              default = TRUE, action ="store", type = "logical", 
              help = "Redo processing of raw data files."),
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
  make_option(c("-z", "--dist_type"), 
              default = "W", action ="store", type = "character", 
              help = "Distribution type from incubation period. One of L, G, W")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

if (!is.null(opt$config)) {
  opt <- yaml::read_yaml(opt$config)
}

# Run checks and print for logs
run_checks_options(opt)
print(opt)

# Other params
p_target <- c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)
period_thresh <- as.Date("2024-09-01")
sd_prior <- 0.3

# Clinical data -----------------------------------------------------------

if (!file.exists(make_clinical_filename()) | opt$redo_data_processing) {
  
  # Basic clinical data
  raw_clinical <- read_raw_data(file_name = "df_kobo_clinical.csv") %>% 
    # Select columns of interest
    select(register_id, age, age_unit, sex, dt_symptom_onset, 
           dt_rash_onset, dt_hospital_admission, dt_interview,
           dt_fever_onset, hospitalized) %>% 
    simplify_colnames() %>% 
    # Add age category
    add_age_num() %>% 
    add_age_cat(col_age = "age_num") %>% 
    tidy_register_id() %>% 
    # Data cleaning
    mutate(
      across(c("register_id"), as.character),
      sex = factor(sex, levels = get_sex_levels()),
      # Make unique admission column based on inpatient/outpatient status
      dt_admission = case_when(hospitalized == "yes" ~ dt_hospital_admission, 
                               T ~ dt_interview),
      # Make simple age cat
      age_cat_simple = forcats::fct_collapse(age_cat, 
                                             "u15" = c("Under 1 year", "1-4 years", "5-14 years"), 
                                             other_level = "15+"),
      period = factor(ifelse(dt_admission < period_thresh, "early", "late")),
      # Simplify hospitalization
      hospitalized = case_when(hospitalized == "unknown" ~ "no",
                               TRUE  ~ hospitalized) %>% 
        factor(levels = c("yes", "no")) 
    )
  
  # Run checks after processing
  run_checks_clinical_data(dat = raw_clinical)
  
  # Save for further use
  saveRDS(raw_clinical, file = make_clinical_filename())
  
} else {
  raw_clinical <- readRDS(make_clinical_filename())
}

#+ fig.height = 4, fig.width = 8
if (opt$do_plots) {
  # Epi curve by age category
  ggplot(raw_clinical, aes(x = dt_admission, fill = age_cat)) +
    geom_bar(stat = "count") +
    theme_bw() +
    ggthemes::scale_fill_few("Dark") +
    labs(x = "date of admission", y = "daily admissions")
}

#+ fig.height = 4, fig.width = 8
if (opt$do_plots) {
  # Case proportion by epi week
  raw_clinical %>% 
    add_epi_week(date_col = "dt_admission") %>% 
    count(epiweek_date, age_cat) %>% 
    group_by(epiweek_date) %>% 
    mutate(frac = n/sum(n)) %>% 
    ggplot(aes(x = epiweek_date, y = frac, fill = age_cat)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    ggthemes::scale_fill_few("Dark") +
    labs(x = "week of admission", y = "proprortion of admissions")
}

# Test data data ----------------------------------------------------------

if (!file.exists(make_test_data_filename()) | opt$redo_data_processing) {
  
  # Load confirmation data
  test_data <- read_raw_data(file_name = "df_kobo_lab.csv") %>% 
    simplify_colnames() %>% 
    select(register_id, contains("genexpert"), dt_record = dt_interview) %>% 
    tidy_register_id() %>% 
    # !! remove duplicates per register_id
    group_by(register_id) %>% 
    arrange(register_id, dt_record) %>% 
    slice(1) %>% 
    ungroup()
  
  # Run checks after processing
  run_checks_test_data(dat = test_data)
  
  # Save
  saveRDS(test_data, file = make_test_data_filename())
  
} else {
  test_data <- readRDS(make_test_data_filename())
}

# Pivot longer by test type
test_data_long <- test_data %>% 
  rename(dt_genexpert_skin_swab_result = dt_genexpert_swab_result) %>% 
  pivot_longer_parts(patterns = c("skin", "throat", "test", "blood"),
                     keys = c("register_id"),
                     colname_parser = genexp_test_colname_parser) %>% 
  mutate(test_result = factor(test_result, levels = c("positive", "negative", "incertaine"))) %>% 
  arrange(register_id, what) %>% 
  filter(!is.na(test_bool)) %>% 
  # !! Keep only first positive test
  group_by(register_id) %>% 
  arrange(desc(test_bool), test_result) %>% 
  slice(1) %>% 
  ungroup()


#+ fig.height = 7, fig.width = 10
if (opt$do_plots) {
  # Epi curve by age category
  raw_clinical %>% 
    left_join(test_data_long) %>% 
    ggplot(aes(x = dt_admission, fill = test_result)) +
    geom_bar(stat = "count") +
    facet_grid(age_cat ~ . ) +
    theme_bw() +
    scale_fill_manual(values = c("red", "blue", "purple")) +
    labs(x = "date of admission", y = "daily admissions")
}

# Contact data  -----------------------------------------------------------

if (!file.exists(make_contacts_filename()) | opt$redo_data_processing) {
  
  # Load contact data
  contact_data <- read_raw_data(file_name = "df_kobo_contacts.csv") %>% 
    select(register_id, contact_type, exposure_location, relationship,
           dt_recent_contact, contact_occurence) %>% 
    tidy_register_id() %>% 
    mutate(exposure_location_cat = case_when(str_detect(exposure_location, "household") ~ "household",
                                             T  ~ "other"))
  
  # Run checks after processing
  run_checks_contact_data(dat = contact_data)
  
  # Save
  saveRDS(contact_data, file = make_contacts_filename())
  
} else {
  contact_data <- readRDS(make_contacts_filename())
}

#  Combine data ---------------------------------------------------------

if (!file.exists(make_combined_filename()) | opt$redo_data_processing) {
  
  # Note that contact_data can have multiple entries per patient as there may be multiple contacts.
  full_data <- raw_clinical %>% 
    left_join(contact_data, by = c("register_id")) %>% 
    # Add unique register_id for each contact
    group_by(register_id) %>% 
    mutate(contact_id = str_c(register_id, row_number(), sep = "-")) %>% 
    ungroup() %>% 
    left_join(test_data_long %>% 
                rename(what_test = what)) %>% 
    mutate(across(contains("dt"), ~ as.Date(.))) %>% 
    # Define exposure based on age
    mutate(exposure_cat = case_when(str_detect(contact_type, "sexual") ~ "sexual",
                                    str_detect(contact_type, "physical") ~ "physical",
                                    T  ~ "other"),
           exposure_cat_simple = case_when(exposure_cat == "sexual" & age_cat_simple == "15+" ~ "sexual",
                                           TRUE ~ "nonsexual") %>% factor(),
           exposure_relation_simple = case_when(relationship %in% c("sex_partner", "spouse") & age_cat_simple == "15+" ~ "sexual_rel",
                                                TRUE ~ "nonsexual_rel") %>% factor()
    )
  
  # Run checks after processing
  run_checks_full_data(dat = full_data)
  
  # Save
  saveRDS(full_data, file = make_combined_filename())
  
} else {
  full_data <- readRDS(make_combined_filename())
}

#+ fig.height = 8, fig.width = 8
if (opt$do_plots) {
  # Symptom onset and contact dates
  full_data %>% 
    mutate(register_id = factor(register_id) %>% forcats::fct_reorder(dt_admission)) %>% 
    pivot_longer(cols = c("dt_recent_contact", 
                          "dt_symptom_onset",
                          "dt_rash_onset",
                          "dt_admission"),
                 names_to = "what",
                 values_to = "date") %>% 
    filter(date > "2024-05-01") %>% 
    ggplot(aes(y = register_id)) +
    geom_line(aes(x = date, group = register_id), lwd = .4, color = "darkgray") +
    geom_point(aes(x = date, pch = what, color = what), size = 1) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank()) +
    ggthemes::scale_color_few("Dark") +
    labs(x = "date", y = "register_id", title = "Participant data")
}

# Reporting delays --------------------------------------------------------

# Compute raw times from symptom onset/rash to admission
raw_reporting_delays <- full_data %>% 
  pivot_longer(cols = c("dt_symptom_onset",
                        "dt_fever_onset",
                        "dt_rash_onset"),
               names_to = "from_what",
               values_to = "from_date") %>% 
  pivot_longer(cols = c("dt_admission"),
               names_to = "to_what",
               values_to = "to_date") %>% 
  mutate(to_date = as.Date(to_date)) %>% 
  mutate(time_diff = difftime(to_date, from_date, units = "days") %>% 
           as.numeric()) 

#+ fig.height = 5, fig.width = 8
if (opt$do_plots) {
  # distributions of times from symptom/rash onset to admission
  raw_reporting_delays %>% 
    # filter(time_diff >= 0) %>% 
    ggplot(aes(x = time_diff, fill = age_cat)) +
    geom_histogram() +
    facet_grid(from_what ~ ., switch = "y", labeller = labeller(test_result = label_both)) +
    ggthemes::scale_fill_few("Dark") +
    theme_bw() +
    labs(x = "time difference [days]",
         title = "Distributions of times from symptom/rash onset to admission")
}


#+ fig.height = 5, fig.width = 8
if (opt$do_plots) {
  # distributions of times from symptom/rash onset to admission
  raw_reporting_delays %>% 
    left_join(raw_clinical %>% select(register_id, dt_admission)) %>% 
    mutate(month = format(dt_admission, "%Y-%m")) %>% 
    ggplot(aes(x = month, y = time_diff)) +
    geom_boxplot(outliers = F) +
    geom_jitter(height = 0, width = .15, alpha = .3) +
    facet_grid(from_what ~ ., switch = "y", labeller = labeller(test_result = label_both)) +
    ggthemes::scale_fill_few("Dark") +
    theme_bw() +
    labs(x = "admission month",
         y = "time from symptom/rash onset to admission",
         title = "Distributions of times from contacts to symptom/rash onset by admission month")
}

# Raw incubation periods --------------------------------------------------
min_exp_date <- as.Date("2024-05-03")

# Compute raw times from reported contact to symptom onset/rash
raw_incubation_periods <- full_data %>% 
  pivot_longer(cols = c("dt_symptom_onset",
                        "dt_fever_onset",
                        "dt_rash_onset"),
               names_to = "to_what",
               values_to = "to_date") %>% 
  mutate(to_date = as.Date(to_date)) %>% 
  pivot_longer(cols = c("dt_recent_contact"),
               names_to = "from_what",
               values_to = "from_date") %>% 
  mutate(time_diff = difftime(to_date, from_date, units = "days") %>% 
           as.numeric()) 

#+ fig.height = 7, fig.width = 8
if (opt$do_plots) {
  # distributions of times from contacts to symptom/rash onset
  raw_incubation_periods %>% 
    filter(!is.na(contact_occurence)) %>% 
    ggplot(aes(x = time_diff, fill = exposure_location_cat)) +
    # geom_density(alpha = .1) + 
    geom_histogram() +
    facet_grid(test_result ~ to_what, switch = "y", scales = "free_y",
               labeller = labeller(test_result = label_both)) +
    ggthemes::scale_fill_few("Dark") +
    theme_bw() +
    labs(x = "time difference [days]",
         title = "Distributions of times from contacts to symptom/rash onset")
}

#+ fig.height = 7, fig.width = 8
if (opt$do_plots) {
  # Distributions of times from contacts to symptom/rash onset by mont
  raw_incubation_periods %>% 
    filter(!is.na(contact_occurence)) %>% 
    mutate(month = format(dt_admission, "%Y-%m")) %>% 
    ggplot(aes(x = month, y = time_diff)) +
    geom_boxplot(outliers = F) +
    geom_jitter(height = 0, width = .15, alpha = .6) +
    facet_grid(test_result ~ to_what, switch = "y", labeller = labeller(test_result = label_both)) +
    ggthemes::scale_fill_few("Dark") +
    theme_bw() +
    labs(x = "admission month",
         y = "time from contact to symptom/rash onset",
         title = "Distributions of times from contacts to symptom/rash onset by admission month")
}

#+ fig.height = 4, fig.width = 6
if (opt$do_plots) {
  # distributions of times from symptom to rash onset
  full_data %>% 
    filter(!is.na(contact_occurence)) %>% 
    mutate(time_diff = difftime(dt_rash_onset, dt_fever_onset, units = "days") %>% 
             as.numeric()) %>% 
    ggplot(aes(x = time_diff, fill = age_cat)) +
    geom_histogram() +
    facet_wrap(~test_result, labeller = labeller(test_result = label_both),
               scales = "free_y") +
    ggthemes::scale_fill_few("Dark") +
    theme_bw() +
    labs(x = "time difference between rash and fever onset [days]",
         title = "Distributions of times from fever to rash onset")
}


#+ fig.height = 4, fig.width = 6
if (opt$do_plots) {
  # distributions of times for individuals that report a single exposure
  raw_incubation_periods %>% 
    filter(!is.na(contact_occurence), time_diff >= 0) %>% 
    ggplot(aes(x = time_diff, fill = period))+
    geom_histogram() +
    facet_grid(contact_occurence ~ to_what, scales = "free") +
    ggthemes::scale_fill_few("Dark") +
    theme_bw() +
    labs(x = "time difference between exposure and symptom onset onset [days]",
         title = "Distributions of times from exposure to symptom onset")
}


# Prepare data for coarseDataTools ----------------------------------------

# Following methods in Lauer et al. 2020
# As we only have most recent contact, we set EL to one month before the first reported case
min_exp_date <- as.Date("2024-04-01")

dic_data <- raw_incubation_periods %>% 
  select(register_id, dt_admission, contact_id, contact_occurence,
         age_cat, to_what, to_date, from_what, from_date, test_result) %>% 
  make_dic_data(min_exp_date = min_exp_date)

#+ fig.height = 8, fig.width = 10
if (opt$do_plots) {
  # dat_sum <- dic_data %>%
  #   mutate(ELnew = EL-ER,
  #          ERnew = ER-ER,
  #          Emid = (ELnew + ERnew)/2,
  #          SLnew = SL-ER,
  #          SRnew = SR-ER,
  #          Smid = (SLnew + SRnew)/2,
  #          UID = forcats::fct_reorder(register_id, to_date))
  # 
  # ggplot(dat_sum, aes(y=factor(UID))) + 
  #   geom_segment(aes(x=ELnew, xend=ERnew, yend=factor(UID)), 
  #                color="#0072B2", size=2, alpha=.25) +
  #   geom_segment(aes(x=Emid, xend=Smid, yend=factor(UID)), size=0.33, color="#999999", alpha = .6) +
  #   geom_point(aes(x=Emid, y=factor(UID)), size=0.5, color="#0072B2") +
  #   geom_point(aes(x=Smid, y=factor(UID)), size=0.5, color="#CC0000") +
  #   #ggtitle("Exposure and symptom onset windows") +
  #   scale_x_continuous("Days from last possible exposure") +
  #   scale_y_discrete("Case") +
  #   theme_bw() +
  #   theme(axis.text.y = element_blank(),
  #         axis.ticks.y= element_blank(),
  #         axis.text.x=element_text(color="black")) +
  #   facet_grid(is.na(contact_occurence) ~ to_what, scales = "free_y", labeller = label_both)
}

saveRDS(dic_data, file = here("generated_data/all_dic_data.rds"))

# Prepare data for stan ---------------------------------------------------

# Set target symptom column
symp_col <- case_when(
  opt$symptom_type == "rash" ~ "dt_rash_onset",
  opt$symptom_type == "fever" ~ "dt_fever_onset",
  opt$symptom_type == "any" ~ "dt_symptom_onset",
  T ~ NA_character_
)

# Prep data
analysis_data <- full_data %>% 
  arrange(register_id, dt_admission, contact_id) %>% 
  # !! keep only data for recent contacts!
  filter(
    !is.na(dt_recent_contact), 
    !is.na(!!rlang::sym(symp_col)),
    !(!!rlang::sym(symp_col) < dt_recent_contact)
  )

# Apply user-defined filters
if (opt$test_subset != "all") {
  analysis_data <- analysis_data %>% 
    filter(test_result == opt$test_subset)
}

# Data for fit with dic.fit
if (opt$model_version == "dic") {
  
  # Make data to DIC format for coarseDataTools
  dic_data <- analysis_data %>% 
    make_dic_data(min_exp_date = min_exp_date,
                  from_col = "dt_recent_contact",
                  to_col = symp_col)
  
  # Run model
  cd_inc_dat <- dic_data %>%
    mutate(type = as.numeric(S_int==0) + as.numeric(E_int==0)) %>%
    select(EL, ER, SL, SR, type, any_of(opt$group_column)) %>%
    as.data.frame() %>%
    filter(EL >= 0, ER >= 0)
  
  saveRDS(list(dic_data = dic_data,
               min_exp_date = min_exp_date,
               cd_inc_dat = cd_inc_dat),
          file = make_inc_period_dic_data_filename(opt = opt))
  
  if (opt$group_column != "none") {
    u_groups <- unique(analysis_data[[opt$group_column]])
    
    # fit using coarseDataTools
    inc_fit_boot <- map(u_groups, function(x) {
      cd_inc_dat %>% 
        filter_at(opt$group_column, ~ . == x) %>% 
        dic.fit(dist = opt$dist_type, n.boots = 100, ptiles = p_target)
      
    })
  } else {
    inc_fit_boot <- cd_inc_dat %>% 
      dic.fit(dist = opt$dist_type,  n.boots = 100, ptiles = p_target)
  }
  
  # Save
  saveRDS(inc_fit_boot, make_inc_period_dic_fit_filename(opt = opt))
  
} else {
  
  # Variables for stan
  
  # Scalars
  N <- nrow(analysis_data)
  M <- analysis_data %>% distinct(register_id) %>% nrow()
  L <- 0
  O <- 60    # window for incubation period (days)
  Q <- 30    # window for reporting
  
  
  # Covariates (none for now)
  X <- array(0, dim = c(M, 0))
  
  # Mappings
  u_ids <- unique(analysis_data$register_id)
  map_to_id <- make_map(vec = analysis_data$register_id, u_ref = u_ids)
  starts <- make_starts(map_to_id)  
  ends <- make_ends(map_to_id)  
  
  # Number of contacts per patient
  n_contacts <- ends - starts  + 1
  
  # Grouping for stratified estimation of incubation period
  if (opt$group_column == "none") {
    J <- 1
    u_groups <- "all"
    map_to_group <- rep(1, N)
  } else {
    u_groups <- unique(analysis_data[[opt$group_column]])
    J <- length(u_groups)
    map_to_group <- make_map(vec = analysis_data[[opt$group_column]],
                             u_ref = u_groups)
  }
  
  # Get log-normal time to lesion resolution from published estimates
  best_pars <- get_lesion_resolution_pars()
  
  if (opt$do_plots) {
    p <- seq(.999, 0.001, by = -.001)
    q <- qlnorm(p, meanlog = best_pars$par[1], sdlog = exp(best_pars$par[2]))
    
    plot(q, 1-p, type = "l")
  }
  
  log_mu_lesions <- best_pars$par[1]
  sd_lesions <- exp(best_pars$par[2])
  
  # Timings
  ref_time <- as.Date("2023-07-01")
  all_admission_times <- sort(diff_time_days(full_data$dt_admission, ref_time))
  t_s <- diff_time_days(analysis_data[[symp_col]], ref_time)
  t_e_r <- diff_time_days(analysis_data$dt_recent_contact, ref_time) 
  t_e_l <- diff_time_days(analysis_data$dt_recent_contact, ref_time)
  all_times <- seq(min(t_s) - O - Q - 1, max(all_admission_times))
  
  # Case analysis data for model version v1
  group_cols <- c("time")
  if (opt$group_column != "none") {
    group_cols <- c(group_cols, opt$group_column)
    across_cols <- u_groups
  } else {
    across_cols <- "n"
  }
  
  analysis_patient_data <- analysis_data %>% 
    group_by(register_id, dt_admission) %>% 
    slice(1) %>% 
    ungroup()
  
  if (opt$group_column != "none") {
    map_patient_to_group <- make_map(vec = analysis_patient_data[[opt$group_column]],
                                     u_ref = u_groups)
  } else {
    map_patient_to_group <- rep(1, nrow(analysis_patient_data))
  }
  
  
  # All patient data for model v1
  patient_dat <- full_data %>% 
    group_by(register_id, dt_admission) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(time = diff_time_days(dt_admission, ref_time))
  
  # Flag for NAs in data
  na_flag <- -9999
  
  # Make observations dataframe
  y_df <- patient_dat %>% 
    group_by_at(group_cols) %>% 
    summarise(n = n()) %>% 
    ungroup()  %>% 
    {
      dat <- .
      if (opt$group_column != "none") {
        dat %>% 
          pivot_wider(values_from = "n",
                      names_from = opt$group_column)
      } else {
        dat
      }
    } %>% 
    complete(time = all_times) %>% 
    arrange(time) %>% 
    # Fill missing data with NA flags
    mutate(
      across(all_of(across_cols), 
             function(x) {
               if (opt$group_column == "period") {
                 
                 case_when(is.na(x) & time < min(all_admission_times) ~ NA, 
                           cur_column() == "early" & time > diff_time_days(period_thresh, ref_time) ~ na_flag,
                           cur_column() == "late" & time <= diff_time_days(period_thresh, ref_time) ~ na_flag,
                           is.na(x) ~ 0,
                           T ~ x)
               } else {
                 case_when(is.na(x) & time < min(all_admission_times) ~ NA, 
                           is.na(x) ~ 0,
                           T ~ x)
               }
             }
      )
    ) 
  
  
  # Define observation indices
  ind_observations <- list()
  ind_no_observations <- list()
  N_observations <- c()
  N_no_observations <- c()
  
  for (i in 1:length(across_cols)) {
    ind_observations <- append(ind_observations, list(which(y_df[, across_cols[i]] > na_flag)))
    ind_no_observations <- append(ind_no_observations, list(which(is.na(y_df[, across_cols[i]]))))
    N_observations <- c(N_observations, length(ind_observations[[i]]))
    N_no_observations <- c(N_no_observations, length(ind_no_observations[[i]]))
  }
  
  # Pad together
  ind_observations_mat <- matrix(0, nrow = length(across_cols), ncol = max(N_observations))
  ind_no_observations_mat <- matrix(0, nrow = length(across_cols), ncol = max(N_no_observations))
  
  for (i in 1:length(across_cols)) {
    ind_observations_mat[i, 1:N_observations[i]] <- ind_observations[[i]]
    ind_no_observations_mat[i, 1:N_no_observations[i]] <- ind_no_observations[[i]]
  }
  
  # Format for stan
  P <- nrow(y_df)
  y <- y_df[, across_cols] %>% as.matrix()
  y[is.na(y)] <- na_flag
  
  # Reporting delay data
  delays <- analysis_data %>% 
    group_by(register_id, dt_admission) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(dt_reporting = diff_time_days(dt_admission, dt_symptom_onset)) %>% 
    filter(!is.na(dt_reporting), dt_reporting > 0)
  
  N_delays <- nrow(delays)
  
  if (opt$group_column == "none") {
    map_delays_group <- rep(1, N_delays)
  } else {
    map_delays_group <- map_dbl(delays[[opt$group_column]], ~ which(u_groups == .))
  }
  
  # Flag multiple contacts
  t_e_l[analysis_data$contact_occurence != "once"] <- t_e_r[analysis_data$contact_occurence != "once"] - O
  # t_e_l[analysis_data$contact_occurence != "once"] <- pmin(diff_time_days("2024-05-01", ref_time), t_e_r[analysis_data$contact_occurence != "once"] - 60)
  censored <- analysis_data$contact_occurence != "once"
  
  # Index of incubation period distribution
  dist_num <- case_when(opt$dist_type == "L" ~ 0,
                        opt$dist_type == "G" ~ 1, 
                        opt$dist_type == "W" ~ 2, 
                        TRUE ~ NA)
  
  # Specific priors on the sd
  sigma_sd_priors <- c(.5, 3, 3)
  sigma_mu_priors <- c(0, 10, 10)
  
  stan_data <- list(
    N = N,
    M = M,
    J = J,
    L = L,
    X = X,
    P = P,
    Q = Q,
    O = O,
    N_observations = N_observations,
    N_no_observations = N_no_observations,
    N_delays = N_delays,
    censored = censored,
    t_s = t_s,
    t_s_int = t_s - min(all_times) + 1,
    t_e_r = t_e_r,
    t_e_l = t_e_l,
    map_to_group = map_to_group,
    starts = starts,
    ends = ends,
    n_contacts = n_contacts,
    log_mu_lesions = log_mu_lesions,
    sigma_lesions = sd_lesions,
    dist = dist_num,
    n_steps_prior = round(t_e_r - t_e_l),
    n_q_gen = 100,
    debug = 0,
    use_data = opt$use_data,
    log_mu_prior = log(9),
    sd_prior = sd_prior,
    y = y,
    dt_reporting = delays$dt_reporting,
    na_flag = na_flag,
    ind_observations = ind_observations_mat,
    ind_no_observations = ind_no_observations_mat,
    n_days_prior = O,
    p_target = p_target,
    n_target_q_gen = length(p_target),
    incid_by_group = opt$group_column != "period",
    x_r = list(),
    x_i = list(),
    map_delays_group = map_delays_group,
    map_patient_to_group = map_patient_to_group,
    sigma_sd_prior = sigma_sd_priors[dist_num + 1],
    sigma_mu_prior = sigma_mu_priors[dist_num + 1],
    n_untrunc = sum((t_s - t_e_r) > 21),
    truncate = 1
  )
  
  saveRDS(list(stan_data = stan_data,
               u_ids = u_ids,
               u_groups = u_groups,
               analysis_data = analysis_data,
               delays = delays,
               y_df = y_df,
               all_times = all_times),
          file = make_inc_period_stan_data_filename(opt = opt))
}
