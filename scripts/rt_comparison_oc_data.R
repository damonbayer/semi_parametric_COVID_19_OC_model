# Description: Run Rt_estim and Epidemia on Orange County data
library(tidyverse)
library(fs)
library(EpiEstim)
source("src/rt_comparison_functions.R")
oc_data <- read_csv("data/oc_data.csv") %>%
  select(time, total_cases = cases, total_tests = tests) %>%
  filter(time <= 42)

dir_create(path("results", "rt_estim"))

gq <- read_csv("data/simulation/oc_like/true_generated_quantities.csv")

# Run Epidemia -----------------------------------------------
data_length <- dim(oc_data)[1]
rates <- c(7 / gq$dur_latent_days, 7 / gq$dur_infectious_days) # pick the generation time

epidemia_weights <- epidemia_hypoexp(data_length, rates) %>% `/`(., sum(.))
delay_weights <- epidemia_gamma(data_length, 1, 7 / gq$dur_latent_days) %>% `/`(., sum(.)) # pick the delay distribution

date <- seq(ymd("2020-07-04"), ymd("2020-07-04") + ddays(data_length), by = "days")

epidemia_data <- data.frame(
  city = "Irvine",
  cases = c(NA, oc_data$total_cases),
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
  seed = 225
)

args$inf <- epiinf(
  gen = epidemia_weights[1:data_length],
  latent = TRUE,
  prior_aux = normal(10, 2)
)
estimnormal_posterior <- do.call(epim, args)

start_date <- min(oc_data$time)
max_date <- max(oc_data$time)

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

write_csv(estimnormal_posterior_rt, path("results", "rt_estim", "rt_comparison_model_id=estimnormal", ext = "csv"))

# Run Rt_estim gamma -----------------------------------------
# run spline for kappa priors
spline <- run_nb_spline(oc_data)

# calculate kappa priors
kappa <- choose_kappa_params(spline)

# calculate quantile for tests
test_quantile <- quantile(oc_data$total_tests)

# first choose rt starting points using epiestim
logrt_start <- get_logrtstart(oc_data, GI_mean = (gq$dur_latent_days + gq$dur_infectious_days) / 7)

# next choose incidence starting points
incid_start <- 1 / 0.066 * oc_data$total_cases

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
estimgamma_posterior <- fit_estimgamma_model(oc_data,
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

write_csv(summary_rtestimgamma, path("results", "rt_estim", "rt_comparison_model_id=estimgamma", ext = "csv"))


# epiestim ----------------------------------------------------------------
mean_time = gq$dur_latent_days + gq$dur_infectious_days
window = 1
GI_mean = mean_time/7
GI_var = 2*(GI_mean/2)^2

ts <- oc_data$time
ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
te <- ts+(window-1)

estimate_R(
    incid = oc_data$total_cases,
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

epiestim_res <- epiestim_weekly[["R"]] %>%
    dplyr::select(t_start,
                  rt_mean = `Mean(R)`,
                  rt_median = `Median(R)`,
                  rt_CI95l = `Quantile.0.025(R)`,
                  rt_CI95u = `Quantile.0.975(R)`) %>%
    mutate(time  = t_start )

epiestim_res <- read_csv("/Users/damon/Documents/semi_parametric_COVID_19_OC_model/results/rt_estim/rt_comparison_model_id=epiestim.csv")
estimagamma <- read_csv("results/rt_estim/rt_comparison_model_id=estimgamma.csv")


write_csv(epiestim_res, path("results", "rt_estim", "rt_comparison_model_id=epiestim", ext = "csv"))