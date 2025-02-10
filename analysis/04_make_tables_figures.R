# This script makes the figures and table for the incubation period analysis


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(posterior)
library(cowplot)
library(optparse)
library(here)

source("analysis/utils.R")

# User-supplied options (if any)
option_list <- list(
  make_option(c("-r", "--redo_processing"), 
              default = FALSE, action ="store", type = "logical", 
              help = "Redo data processing.")
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

# Files
# param_estimates_file <- here("generated_data/all_param_estimates.rds")
cdfs_file <- here("generated_data/all_cdfs.rds")
target_quantiles_file <- here("generated_data/all_target_quantiles.rds")
quantile_traj_file <- here("generated_data/all_quantile_traj.rds")


# Load the results --------------------------------------------------------

# param_estimates <- readRDS(param_estimates_file) %>% 
#   add_factors_to_results() 

inc_period_cdf <- readRDS(cdfs_file) %>% 
  add_factors_to_results() 

inc_period_q_traj <- readRDS(quantile_traj_file) %>% 
  add_factors_to_results() 

inc_period_q <- readRDS(target_quantiles_file) %>% 
  add_factors_to_results() 

# Figure 1 ----------------------------------------------------------------
colors_test_subset <- function() {
  c("purple", "darkorange")
}

# Plot the quantiles for both test subsets and models
filter_for_fig1 <- function(df) {
  df %>% 
    filter(
      model_version == "v0",
      symp == "rash",
      distribution == "W",
      grouping == "overall")
}

pd <- ggstance::position_dodgev(height = .015)
tps <- c(0.05, 0.25, 0.5, 0.75, 0.95)

p_q_v2 <- inc_period_q_traj %>% 
  filter_for_fig1() %>% 
  mutate(comb_draw = str_c(test_subset, .draw)) %>% 
  ggplot(aes(y = ps)) +
  geom_line(aes(x = value, color = test_subset, group = comb_draw), lwd = .3, alpha = 0.075) +
  geom_point(data = inc_period_q %>% 
               filter_for_fig1() %>% 
               filter(ps %in% tps),
             aes(pch = test_subset, x = mean, y = ps, color = test_subset),
             position = pd) +
  geom_errorbarh(data = inc_period_q %>% 
                   filter_for_fig1() %>% 
                   filter(ps %in% tps),
                 aes(xmin = q2.5, xmax = q97.5, y = ps, color = test_subset),
                 position = pd,
                 height = 0, lwd = .4) +
  theme_bw() +
  labs(x = "days from infection", y = "proportion of infections with rashs") +
  scale_fill_manual(values = colors_test_subset())+
  scale_color_manual(values = colors_test_subset()) +
  coord_cartesian(xlim = c(0, 60))

pd <- ggstance::position_dodgev(height = .3)
p_median_v2 <- inc_period_q %>% 
  filter(model_version == "v0",
         symp == "rash",
         distribution == "W",
         ps == 0.5) %>% 
  ggplot(aes(y = group, x = mean, color = test_subset)) +
  geom_vline(data = inc_period_q %>% 
               filter_for_fig1() %>% 
               filter(ps == 0.5) %>% 
               select(-grouping),
             aes(xintercept = mean, color = test_subset), 
             lty = 2, alpha = .4) +
  geom_point(position = pd, aes(pch = test_subset)) +
  geom_errorbarh(aes(xmin = q2.5, xmax = q97.5), height = 0, position = pd, lwd = .3) +
  facet_grid(grouping ~ ., scales = "free", space = "free", switch = "y") +
  theme_bw() +
  labs(y = "strata", x = "median rash incubation period") +
  scale_color_manual(values = colors_test_subset()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(pch = "none")


# Plot the means by grouping
p_fig1_v2 <- plot_grid(
  p_q_v2 +
    labs(color = "Diagnostic test", pch = "Diagnostic test", fill = "Diagnostic test") +
    theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"),
          legend.position = c(.75, .2)),
  p_median_v2 +
    guides(color = "none", fill = "none") +
    theme(plot.margin = unit(c(1, 1, 1, 1), units = "lines"),
          axis.text.y = element_text(size = 6),
          strip.text.y = element_text(size = 6),
          panel.spacing = unit(.1, "lines")),
  nrow = 1,
  rel_widths = c(1.1, 1),
  labels = "auto"
)

ggsave(p_fig1_v2, filename = "figures/inc_period_figure_1_v2.png", 
       width = 10, height = 5.1, dpi = 300)

# Table quantiles ---------------------------------------------------------

target_quantiles <- inc_period_q %>% 
  filter(
    model_version == "v0",
    symp == "rash",
    grouping == "overall",
    ps %in% c(0.05, 0.25, 0.5, 0.75, 0.95)) %>% 
  rename(lo = q2.5, hi = q97.5) %>% 
  make_estimate_text() %>% 
  select(test_subset, distribution_label, ps, txt) %>% 
  mutate(ps = str_c(formatC(ps*100, digits = 0, format = "f"), "%")) %>% 
  pivot_wider(names_from = "ps",
              values_from = "txt") %>% 
  arrange(test_subset, distribution_label)

target_quantiles %>% 
  flextable::as_flextable(show_coltype = FALSE) %>% 
  flextable::fit_to_width(8) %>% 
  flextable::save_as_docx(path = "figures/table_1.docx")

# Supp fig: groupings for fever and any symptom ---------------------------

pd <- ggstance::position_dodgev(height = .3)
p_median_all <- inc_period_q %>% 
  filter(model_version == "v0",
         distribution == "W",
         ps == 0.5) %>% 
  ggplot(aes(y = group, x = mean, color = test_subset)) +
  geom_vline(data = inc_period_q %>% 
               filter(model_version == "v0",
                      distribution == "L",
                      grouping == "overall",
                      ps == 0.5) %>% 
               select(-grouping),
             aes(xintercept = mean, color = test_subset), 
             lty = 2, alpha = .4) +
  geom_point(position = pd, aes(pch = test_subset)) +
  geom_errorbarh(aes(xmin = q2.5, xmax = q97.5), height = 0, position = pd, lwd = .4) +
  facet_grid(grouping ~ symp, scales = "free_y", space = "free_y", switch = "y") +
  theme_bw() +
  labs(y = "strata", x = "median rash incubation period") +
  scale_color_manual(values = colors_test_subset()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.spacing = unit(.1, "lines"),
        legend.position = "top") +
  guides(pch = "none") +
  labs(color = "Diagnostic test", pch = "Diagnostic test", fill = "Diagnostic test") 

ggsave(p_median_all, filename = "figures/sfig_groupgin_all_symp.png", width = 9, height= 8)

