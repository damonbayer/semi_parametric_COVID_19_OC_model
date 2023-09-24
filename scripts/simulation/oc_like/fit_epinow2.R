# fit epinow2
library(EpiNow2)
library(tidyverse)
library(fs)
library(sdprisk)
options(mc.cores = parallelly::availableCores())

target_sim_id <- ifelse(length(commandArgs(trailingOnly = TRUE)) == 0, 1, as.integer(commandArgs(trailingOnly = TRUE[1])))
oc_data <- read_csv("data/oc_data.csv")

# Load data for target_sim_id
simulation_data <-
  read_csv("data/simulation/oc_like/simulated_data.csv") %>%
  select(sim_id = iteration, everything(), -chain) %>%
  pivot_longer(-sim_id, names_to = "name_raw") %>%
  filter(sim_id == target_sim_id) %>%
  mutate(
    name = name_raw %>% str_extract("^.+(?=\\[\\d+\\])") %>% str_remove("data_"),
    time = name_raw %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric()
  ) %>%
  select(sim_id, time, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  left_join(oc_data %>%
              select(time, tests)) %>%
  select(time, total_cases = new_cases, tests) %>%
  rename(total_tests = tests)

dir_create(path("results", "simulation", "oc_like", "epinow2"))

gq <- read_csv("data/simulation/oc_like/true_generated_quantities.csv")

# Run Epinow2 -----------------------------------------------


data_length <- dim(simulation_data)[1]
rates <- c(7 / gq$dur_latent_days, 7 / gq$dur_infectious_days) # pick the generation time


date <- seq(ymd("2020-07-04"), ymd("2020-07-04") + ddays(data_length -1), by = "days")


generation_pmf = c(0, epidemia_weights <- epidemia_hypoexp(data_length, rates) %>% `/`(., sum(.))) # pick generation distribution
delay_pmf = delay_weights <- zero_epidemia_gamma(data_length, 1, 7 / gq$dur_latent_days) %>% `/`(., sum(.)) # pick the delay distribution
mean_time = (gq$dur_latent_days + gq$dur_infectious_days)/7
GI_var = 2*(mean_time/2)^2

mean_delay = gq$dur_latent_days/7
delay_var = 2*(mean_delay/2)^2

epinow2_res <- vector(mode='list', length=200)

for (i in 1:200) {
  # Load data for target_sim_id
  simulation_data <-
    read_csv("data/simulation/oc_like/simulated_data.csv") %>%
    select(sim_id = iteration, everything(), -chain) %>%
    pivot_longer(-sim_id, names_to = "name_raw") %>%
    filter(sim_id == i) %>%
    mutate(
      name = name_raw %>% str_extract("^.+(?=\\[\\d+\\])") %>% str_remove("data_"),
      time = name_raw %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric()
    ) %>%
    select(sim_id, time, name, value) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    left_join(oc_data %>%
                select(time, tests)) %>%
    select(time, total_cases = new_cases, tests) %>%
    rename(total_tests = tests)
    epinow2_data <- data.frame(
      confirm = c(simulation_data$total_cases),
      date = date
    )
    
    
    epinow2_res[[i]] <- estimate_infections(
      epinow2_data,
      generation_time = generation_time_opts(mean = mean_time, sd = sqrt(GI_var), mean_sd = 1.5, sd_sd = 1.5, max = 42, fixed = FALSE),
      delays = delay_opts(list(mean = gq$dur_latent_days/7, sd = sqrt(delay_var), mean_sd = 1.5, sd_sd = 1.5, max = 42), fixed = FALSE),
      truncation = trunc_opts(),
      rt = rt_opts(prior = list(mean = 1.363283, sd = 0.2705766)),
      gp = gp_opts(),
      obs = obs_opts(week_effect = FALSE, scale = list(mean = 0.066, sd = 0.05), phi = c(sqrt(10), 0.6549291)),
      stan = stan_opts(),
      horizon = 0,
      CrIs = c(0.5, 0.8, 0.95),
      filter_leading_zeros = TRUE,
      zero_threshold = Inf,
      id = "estimate_infections",
      verbose = interactive())[["summarised"]]

}

rt_truth <-
  read_csv("data/simulation/oc_like/true_generated_quantities.csv") %>%
  select(starts_with("Rₜ_t")) %>%
  pivot_longer(everything(), names_to = "name_raw") %>%
  mutate(
    name = name_raw %>% str_extract("^.+(?=\\[\\d+\\])"),
    time = name_raw %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric() %>% `-`(1)
  ) %>%
  select(time, Rt = value) %>%
  rename(true_rt = Rt)

# exploring just one 
test = epinow2_res[[1]]
test_intervals = test %>% 
                 filter(variable == "R") %>% 
                 dplyr::select(date, median, lower_80, upper_80) %>% 
                 mutate(time = row_number() - 1) %>% 
                 left_join(rt_truth, by = "time")
rt_metrics(test_intervals, median, upper_80, lower_80)
# grab the 80% credible intervals
epinow2_intervals <- map(epinow2_res, ~.x %>% filter(variable == "R") %>% 
                                          dplyr::select(date, median, lower_80, upper_80) %>% 
                                          mutate(time = row_number() - 1)) 


rt_truth <-
  read_csv("data/simulation/oc_like/true_generated_quantities.csv") %>%
  select(starts_with("Rₜ_t")) %>%
  pivot_longer(everything(), names_to = "name_raw") %>%
  mutate(
    name = name_raw %>% str_extract("^.+(?=\\[\\d+\\])"),
    time = name_raw %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric() %>% `-`(1)
  ) %>%
  select(time, Rt = value) %>%
  rename(true_rt = Rt)



rt_metrics <- function(data, value, upper, lower) {
  metric_one <- data %>%
    mutate(
      dev = abs({{ value }} - true_rt),
      CIW = abs({{ upper }} - {{ lower }}),
      envelope = true_rt >= {{ lower }} & true_rt <= {{ upper }}
    ) %>%
    ungroup() %>%
    filter(!is.na(dev)) %>%
    summarise(
      mean_dev = mean(dev),
      MCIW = mean(CIW),
      mean_env = mean(envelope)
    )
  
  metrics_two <- data %>%
    mutate(
      prev_val = lag({{ value }}),
      prev_rt = lag(true_rt),
      sv = abs({{ value }} - prev_val),
      rt_sv = abs(true_rt - prev_rt)
    ) %>%
    filter(!is.na(sv)) %>%
    ungroup() %>%
    summarise(
      MASV = mean(sv),
      true_MASV = mean(rt_sv)
    )
  
  metrics <- cbind(metric_one, metrics_two)
  
  return(metrics)
}


epinow2_rt_metrics <- map(epinow2_intervals, ~.x %>% 
                             left_join(rt_truth) %>%
                             rt_metrics(
                               median,
                               upper_80,
                               lower_80)) %>%
  bind_rows(.id = "sim_id") %>% 
  mutate(method = "Epinow2")

write_rds(epinow2_intervals, path("results", "simulation", "oc_like", "epinow2_sim_rt_intervals.rds"))


write_csv(summary_rtestimgamma, path("results", "simulation", "oc_like", "epinow2_sim_rt_metrics", ext = "csv"))
