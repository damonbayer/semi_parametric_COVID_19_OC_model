library(tidyverse)
library(fs)
source("src/rt_comparison_functions.R")

rt_truth <-
  read_csv("data/simulation/oc_like/true_generated_quantities.csv") %>%
  select(starts_with("Rₜ_t")) %>%
  pivot_longer(everything(), names_to = "name_raw") %>%
  mutate(
    name = name_raw %>% str_extract("^.+(?=\\[\\d+\\])"),
    time = name_raw %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric() %>% `-`(1)
  ) %>%
  select(time, Rt = value)

full_model_cis <-
  tibble(full_path = dir_ls("results/simulation/oc_like/tidy_posterior_generated_quantities")) %>%
  mutate(sim_id = full_path %>%
           str_extract("(?<=sim_id=)\\d+") %>%
           as.numeric()) %>%
  arrange(sim_id) %>%
  mutate(ci = map(full_path,
                  ~read_csv(.) %>%
                    filter(name == "Rₜ_t",
                           .width == 0.8))) %>% 
  select(-full_path) %>%
  unnest(ci) %>%
  mutate(method = "full")

# ISAAC TO DO: Rt estim and epidemia reuslts ------------------------------
rt_estim_cis <- 
  tibble(full_path = dir_ls("results/simulation/oc_like/estimgamma")) %>%
  mutate(sim_id = full_path %>%
    str_extract("(?<=sim_id=)\\d+") %>%
    as.numeric()) %>%
  arrange(sim_id) %>%
  mutate(ci = map(
    full_path,
    ~ read_csv(.)
  )) %>%
  select(-full_path) %>%
  unnest(ci) %>%
  filter(.width == 0.8) %>%
  mutate(time = time - 1)

epidemia_cis <- tibble(full_path = dir_ls("results/simulation/oc_like/estimnormal")) %>%
  mutate(sim_id = full_path %>%
    str_extract("(?<=sim_id=)\\d+") %>%
    as.numeric()) %>%
  arrange(sim_id) %>%
  mutate(ci = map(
    full_path,
    ~ read_csv(.)
  )) %>%
  select(-full_path) %>%
  unnest(ci) %>%
  filter(.width == 0.8) %>%
  mutate(time = time - 1)

# ISAAC TO DO: compute metrics --------------------------------------------
# Envelope, MCIW, Absolute Deviation, MASV
full_model_metrics <-
  full_model_cis %>%
  left_join(rt_truth, by = "time") %>%
  rename(true_rt = Rt) %>%
  group_by(method, sim_id) %>%
  group_map(~ rt_metrics(.x, .x$value, .x$.upper, .x$.lower)) %>%
  bind_rows(.id = "sim_id") %>%
  mutate(method = "Full Model") %>%
  as_tibble()

estimgamma_metrics <-
  rt_estim_cis %>%
  left_join(rt_truth, by = "time") %>%
  rename(true_rt = Rt) %>%
  group_by(method, sim_id) %>%
  group_map(~ rt_metrics(.x, .x$value, .x$.upper, .x$.lower)) %>%
  bind_rows(.id = "sim_id") %>%
  mutate(method = "Rt-estim-gamma") %>%
  as_tibble()

epidemia_metrics <-
  epidemia_cis %>%
  left_join(rt_truth, by = "time") %>%
  rename(true_rt = Rt) %>%
  group_by(method, sim_id) %>%
  group_map(~ rt_metrics(.x, .x$value, .x$.upper, .x$.lower)) %>%
  bind_rows(.id = "sim_id") %>%
  mutate(method = "Epidemia") %>%
  as_tibble()

all_metrics <-
  bind_rows(full_model_metrics,
            estimgamma_metrics,
            epidemia_metrics) %>%
  mutate(method = fct_inorder(method))

write_csv(all_metrics, "results/simulation/oc_like/simulated_rt_comparison_all_metrics.csv")