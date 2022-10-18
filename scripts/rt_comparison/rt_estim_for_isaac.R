library(tidyverse)
target_sim_id <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 1, as.integer(commandArgs(trailingOnly=T[1])))
source(here::here("src", "damon_functions.R"))
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
  left_join(oc_data %>% select(time, tests))%>% 
  select(time, total_cases = new_cases, tests) %>%
  rename(total_tests = tests)
dir_create(path("results", "simulated_rt_comparison", "estimgamma"))
dir_create(path("results", "simulated_rt_comparison", "estimnormal"))

gq_file_name <- "true_generated_quantities_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv"

gq <- read_csv(here::here("data", "simulated_data", gq_file_name))


# ISAAC TO DO: Run Epidemia -----------------------------------------------

data_length <- dim(simulation_data)[1]
rates <- c(7/gq$dur_latent_days, 7/gq$dur_infectious_days) # pick the generation time


epidemia_weights <- epidemia_hypoexp(data_length, rates)

sum <- sum(epidemia_weights)

epidemia_weights <- epidemia_weights/ sum

delay_weights <- epidemia_gamma(data_length, 1, 7/gq$dur_latent_days) # pick the delay distribution

delay_sum <- sum(delay_weights)

delay_weights <- delay_weights/delay_sum


date <- seq(ymd('2020-07-04'),ymd('2020-07-04') + ddays(data_length), by = "days")

epidemia_data <- data.frame(
  city = "Irvine",
  cases = c(NA, simulation_data$total_cases),
  date = date, 
  day = weekdays(date)
)


rt <- epirt(
  formula = R(city, date) ~ rw(prior_scale = 0.1),
  prior_intercept = normal(log(1.34), 0.2),
  link = 'log'
)

obs <-  epiobs(
  formula = cases ~ 1,
  prior_intercept = rstanarm::normal(location=0.066, scale=0.05),
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
  seed = 225
)


args$inf <- epiinf(gen = epidemia_weights[1:data_length], 
                   latent=TRUE, 
                   prior_aux = normal(10,2))
estimnormal_posterior <- do.call(epim, args)

start_date <- min(simulation_data$time)
max_date <- max(simulation_data$time)


estimnormal_posterior_rt <- data.frame(posterior_rt(estimnormal_posterior)[["draws"]])
names(estimnormal_posterior_rt) <- start_date:(max_date+1)


estimnormal_posterior_rt <- estimnormal_posterior_rt %>%
  mutate(draws = row_number()) %>%
  pivot_longer(!draws, names_to = "epidemia_time", values_to = "value") %>%
  mutate(variable = "rt") %>%
  dplyr::select(variable, epidemia_time, value) %>%
  group_by(variable, epidemia_time) %>%
  median_qi(.width = c(0.5, 0.8, 0.95))

estimnormal_posterior_rt$epidemia_time <- as.numeric(estimnormal_posterior_rt$epidemia_time)

estimnormal_posterior_rt <- estimnormal_posterior_rt%>%
  filter(epidemia_time != 1) %>%
  mutate(time = epidemia_time - 1) %>%
  mutate(method = "estim_normal",
         name = "Rt") %>%
  dplyr::select(time, name, value, .lower, .upper, .width, method)


file_name <- paste0("simulated_rt_comparison_model_id=estimnormal_sim_id=", target_sim_id, ".csv")
write_csv(estimnormal_posterior_rt, here::here("results", "simulated_rt_comparison", "estimnormal", file_name))

# ISAAC TO DO: Run Rt_estim gamma -----------------------------------------
# run spline for kappa priors
spline <- run_nb_spline(simulation_data)


# calculate kappa priors
kappa <- choose_kappa_params(spline)


# calculate quantile for tests
test_quantile <- quantile(simulation_data$total_tests)


# first choose rt starting points using epiestim
logrt_start <- get_logrtstart(simulation_data, GI_mean = (gq$dur_latent_days + gq$dur_infectious_days)/7)

#next choose incidence starting points
incid_start <- 1/0.066 * simulation_data$total_cases

init_func <- function() list(log_incid_rate_raw = 0,
                             log_rt0_raw = 0,
                             rho = 0.066/test_quantile[2],
                             kappa = kappa$par[1],
                             seed_incid_one_raw =1,
                             incid = incid_start,
                             log_rt = logrt_start)
# fit model
estimgamma_posterior <- fit_estimgamma_model(simulation_data, 
                                             gen_params = c(7/gq$dur_latent_days, 7/gq$dur_infectious_days), # pick generation time params
                                             delay_params = c(1, 7/gq$dur_latent_days), # pick delay time params
                                             prev_vals = 4, 
                                             log_rho_mean = log(0.066/test_quantile[2]),
                                             log_rho_sd = 0.3, 
                                             log_r0_mean = log(1.344064),
                                             log_r0_sd = 0.2,
                                             kappa_mean = kappa$par[1],
                                             kappa_sd = kappa$par[2],
                                             iterations = 4000,
                                             thin = 4,
                                             init_func = init_func,
                                             gen_dist = "hypo-exp", # pick the form of generation time
                                             delay_dist = "gamma") # pick the form of delay distribution

# process posterior
summary_rtestimgamma <- summarise_rt_estimgamma(estimgamma_posterior, start_date = 1) %>%
  rename(time = date,
         value = rt) %>%
  mutate(method = "estim_gamma",
         name = "Rt") %>%
  dplyr::select(time, name, value, .lower, .upper, .width, method)

file_name <- paste0("simulated_rt_comparison_model_id=estimgamma_sim_id=", target_sim_id, ".csv")
write_csv(summary_rtestimgamma, here::here("results", "simulated_rt_comparison", "estimgamma", file_name))

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



