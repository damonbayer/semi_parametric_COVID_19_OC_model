# Compute coverage and shrinkage properties for simulated data model fits
library(tidyverse)
library(tidybayes)
library(fs)
source("src/binom.confint.R")

sim_data_dir <- path("data", "simulation", "oc_like")
sim_results_dir <- path("results", "simulation", "oc_like")

sensitivity_model_table <- read_csv("scripts/simulation/oc_like/sensitivity_analysis/sensitivity_model_table.csv")

true_generated_quantities <-
  path(sim_data_dir, "true_generated_quantities.csv") %>%
  read_csv() %>%
  select(-iteration, -chain) %>%
  pivot_longer(everything(), values_to = "true_value") %>%
  separate(col = name,
           into = c("name", "index"),
           sep = "\\[|\\]",
           remove = T,
           fill = "right",
           extra = "drop",
           convert = T) %>%
  mutate(time = case_when(
    name %in% c("α_t", "cases_bb_mean", "cases_nb_mean", "cases_mean", "deaths_mean") ~ index,
    name == "seroprev_mean" ~ 20L,
    TRUE ~ index - 1L)) %>%
  select(name, time, true_value) %>%
  arrange(name, time)

prior_generated_quantities_sd <-
  path(sim_results_dir, "prior_generated_quantities_sd", "prior_generated_quantities_sd.csv") %>%
  read_csv() %>%
  rename(prior_sd = sd)

shrinkage_results <-
  tibble(file_path = dir_ls(path(sim_results_dir, "posterior_generated_quantities_sd"))) %>%
  mutate(sim_id = file_path %>% str_extract("(?<=sim_id=)\\d+") %>% as.integer()) %>%
  arrange(sim_id) %>%
  left_join(sensitivity_model_table) %>% 
  filter(model_design == 1) %>% # all other coverage results would be incorrect because of the prior sd
  mutate(generated_quantities_sd = map(file_path, read_csv)) %>%
  select(-file_path) %>%
  unnest(generated_quantities_sd) %>%
  rename(posterior_sd = sd) %>%
  left_join(prior_generated_quantities_sd) %>%
  mutate(shrinkage = 1 - posterior_sd / prior_sd) %>%
  select(name, time, model_design, shrinkage) %>%
  group_by(name, time, model_design) %>%
  median_qi()


binom_ci_conf_level <- 0.95

coverage_results <-
  tibble(file_path = dir_ls(path(sim_results_dir, "tidy_posterior_generated_quantities"))) %>%
  mutate(sim_id = file_path %>% str_extract("(?<=sim_id=)\\d+") %>% as.integer()) %>%
  arrange(sim_id) %>%
  left_join(sensitivity_model_table) %>% 
  mutate(tidy_posterior_generated_quantities = map(file_path, read_csv))  %>%
  unnest(tidy_posterior_generated_quantities) %>% 
  left_join(true_generated_quantities) %>% 
  filter(!is.na(true_value)) %>% 
  mutate(lower_error = .lower > true_value,
         upper_error = true_value > .upper,
         covered = !lower_error & !upper_error) %>%
  select(-ends_with("_error")) %>%
  group_by(name, time, .width, model_design) %>%
  summarize(x = sum(covered),
            n = n()) %>% 
  mutate(binom.confint(x = x,
                       n = n,
                       methods = "wilson",
                       conf.level = binom_ci_conf_level)) %>%
  mutate(ci_conf_level = binom_ci_conf_level) %>%
  select(name, time, model_design, .width, mean, lower, upper, ci_conf_level)


write_csv(shrinkage_results, path(sim_results_dir, "shrinkage_results.csv"))
write_csv(coverage_results, path(sim_results_dir, "coverage_results.csv"))
