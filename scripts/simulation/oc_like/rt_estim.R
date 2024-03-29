# Run Rt_estim and epidemia on simulated data
library(tidyverse)
library(fs)
target_sim_id <- ifelse(length(commandArgs(trailingOnly = TRUE)) == 0, 1, as.integer(commandArgs(trailingOnly = TRUE[1])))
source("src/rt_comparison_functions.R")
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

dir_create(path("results", "simulation", "oc_like", "estimgamma"))
dir_create(path("results", "simulation", "oc_like", "estimnormal"))

gq <- read_csv("data/simulation/oc_like/true_generated_quantities.csv")

# Run Epidemia -----------------------------------------------
data_length <- dim(simulation_data)[1]
rates <- c(7 / gq$dur_latent_days, 7 / gq$dur_infectious_days) # pick the generation time

epidemia_weights <- epidemia_hypoexp(data_length, rates) %>% `/`(., sum(.))
delay_weights <- epidemia_gamma(data_length, 1, 7 / gq$dur_latent_days) %>% `/`(., sum(.)) # pick the delay distribution

date <- seq(ymd("2020-07-04"), ymd("2020-07-04") + ddays(data_length), by = "days")

epidemia_data <- data.frame(
  city = "Irvine",
  cases = c(NA, simulation_data$total_cases),
  date = date,
  day = weekdays(date)
)

rt <- epirt(
  formula = R(city, date) ~ rw(prior_scale = 0.1),
  prior_intercept = normal(log(1.34), 0.2),
  link = "log"
)

obs <- epiobs(
  formula = cases ~ 1,
  prior_intercept = rstanarm::normal(location = 0.066, scale = 0.05),
  link = "identity",
  i2o = delay_weights[1:data_length]
)
args <- list(
  rt = rt,
  inf = epiinf(gen = epidemia_weights[1:data_length]),
  obs = obs,
  data = epidemia_data,
  iter = 4000,
  thin = 4,
  seed = target_sim_id
)

args$inf <- epiinf(
  gen = epidemia_weights[1:data_length],
  latent = TRUE,
  prior_aux = normal(10, 2)
)
estimnormal_posterior <- do.call(epim, args)

start_date <- min(simulation_data$time)
max_date <- max(simulation_data$time)

estimnormal_posterior_rt <-
  posterior_rt(estimnormal_posterior)[["draws"]] %>%
  data.frame() %>%
  `colnames<-`(start_date:(max_date + 1)) %>%
  mutate(draws = row_number()) %>%
  pivot_longer(!draws,
    names_to = "epidemia_time",
    values_to = "value",
    names_transform = list(epidemia_time = as.integer)
  ) %>%
  mutate(variable = "rt") %>%
  dplyr::select(variable, epidemia_time, value) %>%
  group_by(variable, epidemia_time) %>%
  median_qi(.width = c(0.5, 0.8, 0.95)) %>%
  filter(epidemia_time != 1) %>%
  mutate(time = epidemia_time - 1) %>%
  mutate(
    method = "estim_normal",
    name = "Rt"
  ) %>%
  dplyr::select(time, name, value, .lower, .upper, .width, method)

estimnormal_file_name <- paste0("simulated_rt_comparison_model_id=estimnormal_sim_id=", target_sim_id, ".csv")
write_csv(estimnormal_posterior_rt, path("results", "simulation", "oc_like", "estimnormal", estimnormal_file_name))

# Run Rt_estim gamma -----------------------------------------
# run spline for kappa priors
spline <- run_nb_spline(simulation_data)

# calculate kappa priors
kappa <- choose_kappa_params(spline)

# calculate quantile for tests
test_quantile <- quantile(simulation_data$total_tests)

# first choose rt starting points using epiestim
logrt_start <- get_logrtstart(simulation_data, GI_mean = (gq$dur_latent_days + gq$dur_infectious_days) / 7)

# next choose incidence starting points
incid_start <- 1 / 0.066 * simulation_data$total_cases

init_func <- function() {
  list(
    log_incid_rate_raw = 0,
    log_rt0_raw = 0,
    rho = 0.066 / test_quantile[2],
    kappa = kappa$par[1],
    seed_incid_one_raw = 1,
    incid = incid_start,
    log_rt = logrt_start
  )
}

# fit model
estimgamma_posterior <- fit_estimgamma_model(simulation_data,
  gen_params = c(7 / gq$dur_latent_days, 7 / gq$dur_infectious_days), # pick generation time params
  delay_params = c(1, 7 / gq$dur_latent_days), # pick delay time params
  prev_vals = 4,
  log_rho_mean = log(0.066 / test_quantile[2]),
  log_rho_sd = 0.3,
  log_r0_mean = log(1.344064),
  log_r0_sd = 0.2,
  kappa_mean = kappa$par[1],
  kappa_sd = kappa$par[2],
  iterations = 4000,
  thin = 4,
  init_func = init_func,
  gen_dist = "hypo-exp", # pick the form of generation time
  delay_dist = "gamma"
) # pick the form of delay distribution

# process posterior
summary_rtestimgamma <-
  summarise_rt_estimgamma(estimgamma_posterior, start_date = 1) %>%
  rename(
    time = date,
    value = rt
  ) %>%
  mutate(
    method = "estim_gamma",
    name = "Rt"
  ) %>%
  dplyr::select(time, name, value, .lower, .upper, .width, method)

estimgamma_file_name <- paste0("simulated_rt_comparison_model_id=estimgamma_sim_id=", target_sim_id, ".csv")
write_csv(summary_rtestimgamma, path("results", "simulation", "oc_like", "estimgamma", estimgamma_file_name))
