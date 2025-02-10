# This script makes the configs to run all analyses

# Preamble ----------------------------------------------------------------
library(here)
library(purrr)
library(tibble)
library(stringr)
library(dplyr)

# Write configs -----------------------------------------------------------

for (subset in c("overall", "groupings")) {
  
  if (subset == "overall") {
    versions <- c("v0", "v1")
  } else {
    versions <- c("v0")
  }
  
  for (version in versions) {
    
    dir.create(str_glue("analysis/configs/{subset}_{version}"))
    
    # Make all combinations of specs
    config_specs <- expand.grid(
      model_version = version,
      test_subset = c("positive", "all"),
      dist_type = c("L", "G", "W"),
      symptom_type = c("rash", "fever", "any"),
      group_column = c("none", "age_cat_simple", "sex", "hospitalized", "period", "exposure_cat_simple", "exposure_relation_simple")
    ) %>% 
      as_tibble()
    
    if (subset == "overall") {
      config_specs <- config_specs %>% 
        filter(group_column == "none")
    } else {
      config_specs <- config_specs %>% 
        filter(group_column != "none")
    }
    
    # Write configs
    walk(1:nrow(config_specs), function(x){
      config <- yaml::read_yaml(here("analysis/config_defaults.yml"))
      
      this_specs <- unlist(config_specs[x, ])
      
      for (i in names(this_specs)) {
        config[[i]] <- this_specs[i]
      }
      
      yaml::write_yaml(config, here(str_glue("analysis/configs/{subset}_{version}/config_{subset}_{version}_{x}.yml")))
      
    })
  }
}

