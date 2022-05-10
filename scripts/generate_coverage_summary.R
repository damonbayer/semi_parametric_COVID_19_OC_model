library(tidyverse)
library(tidybayes)
library(fs)

true_parameters <-
  read_csv("data/simulated_data/true_parameters_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  select(-iteration, -chain) %>% 
  pivot_longer(everything(), values_to = "true_value")


coverage_summary <- 
  tibble(full_path = dir_ls("results/simulated_posterior_samples_csv")) %>% 
  mutate(tidy_summary = map(full_path,
                            ~read_csv(.) %>% 
                              select(-iteration, -chain) %>% 
                              pivot_longer(everything()) %>% 
                              group_by(name) %>% 
                              median_qi(.width = c(0.5, 0.8, 0.95)))) %>% 
  unnest(tidy_summary) %>% 
  select(-full_path) %>% 
  left_join(true_parameters) %>% 
  mutate(lower_error = .lower > true_value,
         upper_error = true_value > .upper) %>% 
  mutate(covered = !lower_error & !upper_error) %>% 
  select(name, .width, lower_error, upper_error, covered) %>% 
  pivot_longer(cols = c(lower_error, upper_error, covered), names_to = "type") %>% 
  group_by(name, .width, type) %>% 
  summarize(value = mean(value)) %>% 
  pivot_wider(names_from = type, values_from = value)

write_csv(coverage_summary, "results/coverage_summary.csv")
