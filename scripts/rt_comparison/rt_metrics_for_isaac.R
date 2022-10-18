library(tidyverse)
library(fs)

rt_truth <- 
  read_csv("data/simulated_data/true_generated_quantities_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  select(starts_with("Rₜ_t")) %>% 
  pivot_longer(everything(), names_to = "name_raw") %>% 
  mutate(name = name_raw %>% str_extract("^.+(?=\\[\\d+\\])"),
         time = name_raw %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  select(time, Rt = value)


# full_model_cis <- 
#   tibble(full_path = dir_ls("results/simulated_generated_quantities_summary")) %>% 
#   filter(full_path %>% str_detect("prior", negate = T)) %>% 
#   mutate(sim_id = full_path %>% 
#            str_extract("(?<=seed=)\\d+") %>% 
#            as.numeric()) %>% 
#   arrange(sim_id) %>% 
#   mutate(ci = map(full_path,
#                   ~read_csv(.) %>% 
#                     select(name_raw = parameters,
#                            .lower = `10.0%`,
#                            .upper = `90.0%`) %>% 
#                     filter(name_raw %>% str_starts("Rₜ_t")) %>% 
#                     mutate(name = name_raw %>% str_extract("^.+(?=\\[\\d+\\])") %>% str_remove("data_"),
#                            time = name_raw %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
#                     select(time, .lower, .upper))) %>% 
#   select(-full_path) %>% 
#   unnest(ci) %>% 
#   mutate(.width = 0.8,
#          method = "full")



# ISAAC TO DO: Rt estim and epidemia reuslts ------------------------------
rt_estim_cis <- tibble(full_path = dir_ls("results/simulated_rt_comparison/estimgamma")) %>% 
  filter(full_path %>% str_detect("prior", negate = T)) %>%
  mutate(sim_id = full_path %>%
           str_extract("(?<=sim_id=)\\d+") %>%
           as.numeric()) %>%
  arrange(sim_id) %>%
  mutate(ci = map(full_path,
                  ~read_csv(.))) %>%
  select(-full_path) %>%
  unnest(ci) %>%
  filter(.width == 0.8) %>%
  mutate(time = time - 1)

epidemia_cis <- tibble(full_path = dir_ls("results/simulated_rt_comparison/estimnormal")) %>% 
  filter(full_path %>% str_detect("prior", negate = T)) %>%
  mutate(sim_id = full_path %>%
           str_extract("(?<=sim_id=)\\d+") %>%
           as.numeric()) %>%
  arrange(sim_id) %>%
  mutate(ci = map(full_path,
                  ~read_csv(.))) %>%
  select(-full_path) %>%
  unnest(ci) %>%
  filter(.width == 0.8) %>% 
  mutate(time = time - 1)

# ISAAC TO DO: compute metrics --------------------------------------------
# Envelope, MCIW, Absolute Deviation, MASV
full_model_metrics <- full_model_cis %>%
  left_join(rt_truth, by = "time") %>%
  rename(true_rt = Rt) %>%
  group_by(method, sim_id) %>%
  group_map(~rt_metrics(.x, .x$value, .x$.upper, .x$.lower)) %>%
  bind_rows(.id = "sim_id") %>%
  mutate(method = "full_model")

estimgamma_metrics <- rt_estim_cis %>% 
  left_join(rt_truth, by = "time") %>%
  rename(true_rt = Rt) %>%
  group_by(method, sim_id) %>%
  group_map(~rt_metrics(.x, .x$value, .x$.upper, .x$.lower)) %>%
  bind_rows(.id = "sim_id") %>%
  mutate(method = "estim_gamma")

epidemia_metrics <- epidemia_cis %>% 
  left_join(rt_truth, by = "time") %>%
  rename(true_rt = Rt) %>%
  group_by(method, sim_id) %>%
  group_map(~rt_metrics(.x, .x$value, .x$.upper, .x$.lower)) %>%
  bind_rows(.id = "sim_id") %>%
  mutate(method = "epidemia")

write_csv(full_model_metrics, here::here("results", "full_model_metrics.csv"))
write_csv(estimgamma_metrics, here::here("results", "estimgamma_metrics.csv"))
write_csv(epidemia_metrics, here::here("results", "epidemia_metrics.csv"))
