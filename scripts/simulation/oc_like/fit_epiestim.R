# fit epiestim to simulated data
library(tidyverse)
library(fs)
library(EpiEstim)
target_sim_id <- ifelse(length(commandArgs(trailingOnly = TRUE)) == 0, 1, as.integer(commandArgs(trailingOnly = TRUE[1])))
# source("src/rt_comparison_functions.R")
oc_data <- read_csv("data/oc_data.csv")


gq <- read_csv("data/simulation/oc_like/true_generated_quantities.csv")

# Run epiestim -----------------------------------------------
mean_time = gq$dur_latent_days + gq$dur_infectious_days
epiestim_res <- vector(mode='list', length=200)

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
  # fit epiestim
  # fit model  
  window = 1
  GI_mean = mean_time/7
  GI_var = 2*(GI_mean/2)^2
  
  ts <- simulation_data$time
  ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
  te <- ts+(window-1)
  
  
  estimate_R(
    incid = simulation_data$total_cases,
    method = "uncertain_si",
    config = make_config(
      list(
        mean_si = GI_mean,
        min_mean_si = 1,
        max_mean_si = GI_mean + 1,
        std_mean_si = 1.5,
        std_std_si = 1.5,
        std_si = sqrt(GI_var),
        min_std_si = sqrt(GI_var)*.8,
        max_std_si = sqrt(GI_var)*1.2,
        n1 = 50,
        n2 = 100, 
        t_start=ts,
        t_end=te
      )
    )
  ) -> epiestim_weekly
  
  epiestim_res[[i]] <- epiestim_weekly[["R"]] %>%
    dplyr::select(t_start, 
                  rt_mean = `Mean(R)`, 
                  rt_median = `Median(R)`,
                  rt_CI95l = `Quantile.0.025(R)`,
                  rt_CI95u = `Quantile.0.975(R)`) %>%
    mutate(time  = t_start -1) 
  
  
} 


# processing metrics ------------------------------------------------------


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

epiestim_rt_metrics <- map(epiestim_res, ~.x %>% 
                                          left_join(rt_truth) %>%
                                  rt_metrics(
                                  rt_median,
                                  rt_CI95u,
                                  rt_CI95l)) %>%
  bind_rows(.id = "sim_id") %>% 
  mutate(method = "EpiEstim")

write_csv(epiestim_rt_metrics, here::here("scripts", "simulation", "oc_like", "epiestim_sim_rt_metrics.csv"))

