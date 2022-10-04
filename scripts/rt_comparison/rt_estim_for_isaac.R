library(tidyverse)
target_sim_id <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 1, as.integer(commandArgs(trailingOnly=T[1])))

oc_data <- read_csv("data/oc_data.csv")

# Load data for target_sim_id 
simulation_data <-
  read_csv("data/simulated_data/simulated_data_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>%
  select(sim_id = iteration, everything(), -chain) %>%
  pivot_longer(-sim_id, names_to = "name_raw") %>% 
  filter(sim_id == target_sim_id) %>% 
  mutate(name = name_raw %>% str_extract("^.+(?=\\[\\d+\\])") %>% str_remove("data_"),
         time = name_raw %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  select(sim_id, time, name, value) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  left_join(oc_data %>% select(time, tests)) %>% 
  select(time, cases = new_cases, tests)



# ISAAC TO DO: Run Epidemia -----------------------------------------------

# ISAAC TO DO: Run Rt_estim gamma -----------------------------------------

# ISAAC TO DO: Save 80% CI's somewhere ------------------------------------
# Yes, 80% CI, not 95% CI

# I would probably save them in results/simulated_rt_comparison
# with file names simulated_rt_comparison_sim_id=xxx.csv
# I suggest this format, but you can do whatever you want.
# Not sure if you estimate Rt at 0 or not.

#  time name  value .lower .upper .width method     
# <dbl> <chr> <dbl>  <dbl>  <dbl>  <dbl> <chr>      
#     0 Rₜ_t   1.27   1.18   1.36    0.5 method_name
#     1 Rₜ_t   1.22   1.17   1.29    0.5 method_name
#     2 Rₜ_t   1.34   1.29   1.40    0.5 method_name
#     3 Rₜ_t   1.40   1.34   1.51    0.5 method_name
#     4 Rₜ_t   1.34   1.28   1.41    0.5 method_name
#     5 Rₜ_t   1.31   1.24   1.36    0.5 method_name
#     6 Rₜ_t   1.34   1.30   1.38    0.5 method_name
#     7 Rₜ_t   1.26   1.21   1.30    0.5 method_name
#     8 Rₜ_t   1.24   1.16   1.33    0.5 method_name
#     9 Rₜ_t   1.24   1.19   1.29    0.5 method_name




# Compare Epidemia, and Rt_estim gamma and my method

# Step 1.
# Generate 95% confidence/credible intervals for 200 simulated data sets.
# That gets stored somewhere in a format that Damon understands.
# Isaac can just chack for 10 models or whatever locally.

# Step 2.
# Damon will add semi-parametric 95% credible intervals to that table

# Step 3.
# Isaac has code that turns that table into useful metrics.
# Envelope, MCIW, Absolute Deviation, MASV



current_sim_id <- 1
# For each sim_id (map or loop or whatever):
# current_dat <- 
  dat %>%
  filter(sim_id == current_sim_id) %>% 
    select(week = time, cases = new_cases, tests = tests) %>% 
    dput()

current_dat %>% 
  select(-sim_id, -seroprev_cases) %>% 
  rename_all(~str_remove(., "new_"))



