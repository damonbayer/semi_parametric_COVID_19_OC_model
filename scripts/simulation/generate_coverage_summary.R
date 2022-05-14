library(tidyverse)
library(tidybayes)
library(fs)

true_parameters <-
  read_csv("data/simulated_data/true_parameters_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  select(-iteration, -chain) %>% 
  pivot_longer(everything(), values_to = "true_value")

true_gq <-
  read_csv("data/simulated_data/true_generated_quantities_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  select(-iteration, -chain) %>% 
  pivot_longer(everything(), values_to = "true_value")

prior_parameters_std <- 
  read_csv("results/simulated_posterior_samples_summary/simulated_prior_samples_summary_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  select(name = parameters, prior_std = std)

prior_gq_std <- 
  read_csv("results/simulated_generated_quantities_summary/simulated_prior_generated_quantities_summary_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  select(name = parameters, prior_std = std)


posterior_samples_coverage_summary <- 
  tibble(full_path = dir_ls("results/simulated_posterior_samples_summary")) %>% 
  filter(str_detect(full_path, "prior", negate = T)) %>% 
  mutate(summary = map(full_path, read_csv)) %>% 
  mutate(sim_id = full_path %>% 
           str_extract("(?<=seed=)\\d+") %>% 
           as.integer()) %>% 
  select(-full_path) %>%
  arrange(sim_id) %>% 
  unnest(summary) %>% 
  rename(name = parameters,
         .lower = `10.0%`,
         .upper = `90.0%`) %>% 
  left_join(true_parameters) %>% 
  left_join(prior_parameters_std) %>% 
  mutate(lower_error = .lower > true_value,
         upper_error = true_value > .upper,
         covered = !lower_error & !upper_error,
         shrinkage = 1 - std / prior_std) %>% 
  select(sim_id, name, lower_error, upper_error, covered, shrinkage) %>% 
  pivot_longer(-c(sim_id, name), names_to = "stat") %>% 
  group_by(name, stat) %>% 
  summarize(value = mean(value), .groups = "drop") %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  mutate(time = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  mutate(name = if_else(str_detect(name,"^.+\\["), str_extract(name, "^.+(?=\\[)"), name)) %>% 
  select(name, time, shrinkage, covered, lower_error, upper_error)


generated_quantities_coverage_summary <- 
  tibble(full_path = dir_ls("results/simulated_generated_quantities_summary/")) %>% 
  filter(str_detect(full_path, "prior", negate = T)) %>% 
  mutate(summary = map(full_path, read_csv)) %>% 
  mutate(sim_id = full_path %>% 
           str_extract("(?<=seed=)\\d+") %>% 
           as.integer()) %>% 
  select(-full_path) %>%
  arrange(sim_id) %>% 
  unnest(summary) %>% 
  rename(name = parameters,
         .lower = `10.0%`,
         .upper = `90.0%`) %>% 
  left_join(true_gq) %>% 
  left_join(prior_gq_std) %>% 
  mutate(lower_error = .lower > true_value,
         upper_error = true_value > .upper,
         covered = !lower_error & !upper_error,
         shrinkage = 1 - std / prior_std) %>% 
  select(sim_id, name, lower_error, upper_error, covered, shrinkage) %>% 
  pivot_longer(-c(sim_id, name), names_to = "stat") %>% 
  group_by(name, stat) %>% 
  summarize(value = mean(value), .groups = "drop") %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  mutate(time = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  mutate(name = if_else(str_detect(name,"^.+\\["), str_extract(name, "^.+(?=\\[)"), name)) %>% 
  select(name, time, shrinkage, covered, lower_error, upper_error)

write_csv(posterior_samples_coverage_summary, "results/posterior_samples_coverage_summary.csv")
write_csv(generated_quantities_coverage_summary, "results/generated_quantities_coverage_summary.csv")
