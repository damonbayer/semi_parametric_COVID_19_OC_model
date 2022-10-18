# damon functions
library(tidyverse)
library(EpiEstim)
library(tidybayes)
library(truncnorm)
library(rstan)
library(lubridate)
library(brms)
library(epidemia)
library(rstanarm)
library(sdprisk)

# Create automated process for choosing kappa parameter from a spl --------

run_nb_spline <- function(data, 
                          seed = 17,
                          iter = 4000,
                          warmup = 1000,
                          thin = 10, 
                          refresh = 0,
                          adapt_delta = 0.99) {
  spline_model <- brm(bf(total_cases ~ s(time)),
                      data = data, family = negbinomial(), cores = 4, seed = seed,
                      iter = iter, warmup = warmup, thin = thin, refresh = refresh,
                      control = list(adapt_delta = adapt_delta))
  
  return(spline_model)
}

compare_kappa_quantiles <- function(candidate_params, true_quantiles) {
  candidate_quantiles <- quantile(rtruncnorm(10000, 
                                             a = 0, 
                                             mean = candidate_params[1], 
                                             sd = candidate_params[2]),
                                  c(0.025, 0.975))
  loss <- (true_quantiles[1] - candidate_quantiles[1])^2 + (true_quantiles[2] - candidate_quantiles[2])^2
  return(loss)
}


choose_kappa_params <- function(spline_posterior) {
  posterior_pars <- summary(spline_posterior)
  start_mean <- posterior_pars[["spec_pars"]][[1]]
  start_sd <- posterior_pars[["spec_pars"]][[2]]
  true_lb <- posterior_pars[["spec_pars"]][[3]]
  true_ub <- posterior_pars[["spec_pars"]][[4]]
  
  start_params <- c(start_mean, start_sd)
  true_quantiles <- c(true_lb, true_ub)
  
  optim_params <- optim(par = start_params,
                        fn = compare_kappa_quantiles, 
                        true_quantiles = true_quantiles)
}

# fit exp_seed model ------------------------------------------------------
# this is the current version of the model
# assuming generation time is a hypo-expo distribution
# assuming delay time is a gamma distribution, admitting the zero case
fit_estimgamma_model <- function(data,
                                 gen_params,
                                 delay_params,
                                 prev_vals,
                                 log_nu_mean = -2,
                                 log_nu_sd = 0.7,
                                 log_sigma_mean = -0.6,
                                 log_sigma_sd = 0.6,
                                 log_rho_mean,
                                 log_rho_sd,
                                 log_r0_mean = log(1),
                                 log_r0_sd = 0.75,
                                 kappa_mean,
                                 kappa_sd,
                                 init_func,
                                 iterations = 2000,
                                 thin = 2,
                                 adapt_delta = 0.99,
                                 treedepth = 12,
                                 seed = 45,
                                 chain = 4,
                                 gen_dist = "hypo-exp",
                                 delay_dist = "gamma") {
  
  data_length <- dim(data)[1]
  
  if (gen_dist == "hypo-exp") {
    gen_weights <- epidemia_hypoexp(data_length, gen_params)
    
  }
  
  
  
  if (gen_dist == "log-normal") {
    gen_weights <- epidemia_lognormal(data_length, gen_params)
  }
  
  if (gen_dist == "weibull") {
    gen_weights <- epidemia_weibull(data_length, gen_params)
  }
  
  
  if (delay_dist == "gamma") {
    delay_weights <- zero_epidemia_gamma(data_length, 
                                         delay_params[1], 
                                         delay_params[2])
  }
  
  model_object <- list(n = data_length, 
                       d = data_length,
                       w = gen_weights,
                       delay_weights = delay_weights,
                       obs = data$total_cases,
                       test = data$total_tests,
                       prev_vals = 4,
                       log_incid_rate_mean = log_nu_mean,
                       log_incid_rate_sd = log_nu_sd,
                       log_sigma_mu = log_sigma_mean,
                       log_sigma_sd = log_sigma_sd,
                       log_rho_mu = log_rho_mean,
                       log_rho_sd = log_rho_sd,
                       log_r0_mu = log_r0_mean,
                       log_r0_sd = log_r0_sd,
                       kappa_mu = kappa_mean,
                       kappa_sd = kappa_sd)
  
  
  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = treedepth)
  
  model_fit <- stan(file = here::here("src", "rt_estim_gamma.stan"),
                    data = model_object,
                    seed = seed,
                    iter = iterations,
                    thin = thin,
                    chain = chain,
                    init = init_func,
                    control = control_list)
  
  return(model_fit)
}

# create estimgamma rt posterior for sim data, no truth involved --------------------
summarise_rt_estimgamma <- function(stan_posterior,
                                             start_date,
                                             include_chains = c(1,2,3,4)){
  
  
  rt_posterior <- stan_posterior %>%
    spread_draws(log_rt[i]) %>%
    filter(.chain %in% include_chains) %>%
    group_by(.draw) %>% 
    arrange(.draw, i)  %>%
    mutate(date = i + 0 + start_date -1) %>%
    mutate(rt = exp(log_rt)) %>%
    dplyr::select(date, rt) %>%
    group_by(date) %>% 
    median_qi(.width = c(0.5, 0.8, 0.95)) 
  return(rt_posterior)
}


# rt_metrics --------------------------------------------------------------
# operating characteristics

rt_metrics<- function(data, value, upper, lower) {
  metric_one <- data %>%
    mutate(dev = abs({{ value }} - true_rt),
           CIW = abs({{ upper }} - {{ lower }}),
           envelope = true_rt >= {{ lower }} & true_rt <=  {{ upper }}) %>%
    ungroup() %>%
    filter(!is.na(dev)) %>%
    summarise(mean_dev = mean(dev),
              MCIW = mean(CIW),
              mean_env = mean(envelope))
  
  metrics_two <- data %>%
    mutate(prev_val = lag({{ value }}),
           prev_rt = lag(true_rt),
           sv = abs({{ value }} - prev_val),
           rt_sv = abs(true_rt - prev_rt)) %>%
    filter(!is.na(sv)) %>%
    ungroup() %>%
    summarise(MASV = mean(sv),
              true_MASV = mean(rt_sv))
  
  metrics <- cbind(metric_one, metrics_two)
  
  return(metrics)
}


# use epiestim to choose initial conditions for rt ------------------------

get_logrtstart <- function(data,
                           window = 1, 
                           GI_mean = 11.5/7
) {
  
  window = window
  GI_mean = GI_mean
  GI_var = 2*(GI_mean/2)^2
  
  ts <- data$time
  ts <- ts[ts > 1 & ts <= (max(ts)-window+1)]
  te <- ts+(window-1)
  
  estimate_R(
    incid = data$total_cases,
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
  ) -> ee_outs
  
  ee_quantile <- ee_outs[["R"]] %>%
    dplyr::select(t_start, 
                  rt_mean = `Mean(R)`, 
                  rt_median = `Median(R)`,
                  rt_CI95l = `Quantile.0.025(R)`,
                  rt_CI95u = `Quantile.0.975(R)`) %>%
    mutate(time  = t_start) 
  
  
  log_ee_median <- log(ee_quantile %>% pull(rt_median))
  
  first_one <- log_ee_median[1]
  rt_start <- c(first_one, log_ee_median)
  
  return(rt_start)
}

# fit estim_normal model ------------------------------------------------------
# this is the current version of the model
# assuming generation time is a hypo-expo distribution
# assuming delay time is a gamma distribution, no zero case
# gen_params = c(log(7.872346) + log(1/7),
#                0.642713)
# delay_params = c(4.05, 7*0.74)
# iterations = 50
# seed = 12345
# init = FALSE
# thin = 3
# gen_distribution = "log-normal"
# rt_int_mean = 0
# rt_int_scale = 0.75
# obs_int_loc = 0.06
# obs_int_scale = 0.07
# psi_mean = 10
# psi_sd = 2
# fit_estimnormal_model <- function(data,
#                                   gen_params,
#                                   delay_params,
#                                   rt_int_mean = 0,
#                                   rt_int_scale = 0.75,
#                                   obs_int_loc = 0.06,
#                                   obs_int_scale = 0.07,
#                                   psi_mean = 10,
#                                   psi_sd = 2,
#                                   iterations = 2000,
#                                   seed = 45,
#                                   init= FALSE,
#                                   init_func = NULL,
#                                   thin = 1,
#                                   gen_distribution = "hypo-expo",
#                                   delay_distribution = "gamma") {
#   
# 
#   data_length <- dim(data)[1]
#   
#   if (gen_distribution == "hypo-expo") {
#     hypoexpo_weights <- epidemia_hypoexp(data_length, gen_params)
#     
#     sum <- sum(hypoexpo_weights)
#     epidemia_weights <- hypoexpo_weights/sum
#     
#   }
#   
#   if (gen_distribution == "log-normal") {
#     lognormal_weights <- epidemia_lognormal(data_length, gen_params)
#     
#     sum <- sum(lognormal_weights)
#     epidemia_weights <- lognormal_weights/sum
#     
#   }
#   
#   if (delay_distribution == "gamma") {
#     delay_weights <- epidemia_gamma(data_length, delay_params[1], delay_params[2])
#     
#     
#     sum_delay <- sum(delay_weights)
#     
#     delay_weights <- delay_weights/sum_delay
#     
#   }
#   
#   
#   first_date <- as.Date("2020-08-01")
#   date <- rep(first_date, data_length+1)
#   for (i in 2:(data_length+1)){
#     date[i] <- date[i-1] + ddays(1)
#   }
# 
#   data <- data.frame(
#     city = "Everywhere",
#     cases = c(NA, data$total_cases),
#     date = date,
#     day = weekdays(date)
#   )
#   
#   
#   rt <- epirt(
#     formula = R(city, date) ~ rw(prior_scale = 0.1),
#     prior_intercept = normal(log(1), 0.2),
#     link = 'log'
#   )
#   
#   obs <-  epiobs(
#     formula = cases ~ 1,
#     prior_intercept = rstanarm::normal(location=0.02, scale=0.05),
#     link = "identity",
#     i2o = delay_weights[1:data_length]
#     #prior_aux = rstanarm::normal(location = 1/58, scale = 1/35)
#   )
#   
#   
#   if (init == TRUE){
#     args <- list(
#       rt = rt,
#       inf = epiinf(gen = epidemia_weights[1:data_length]),
#       obs = obs,
#       data = data,
#       iter = iterations,
#       seed = seed,
#       init = init_func,
#       thin = thin
#     )
#     
#   }
#   
#   if (init == FALSE){
#     args <- list(
#       rt = rt,
#       inf = epiinf(gen = epidemia_weights[1:data_length]),
#       obs = obs,
#       data = data,
#       iter = iterations,
#       seed = seed,
#       thin = thin
#     )
#     
#   }
#   
# 
#   #to redo with incidence as parameter
#   args$inf <- epiinf(gen = epidemia_weights[1:data_length],
#                      latent=TRUE,
#                      prior_aux = normal(psi_mean,psi_sd))
#   fm2 <- do.call(epim, args)
#   
#   
#   
#   return(fm2)
# }

# Discretize  Distributions ------------------------------------------
# Epidemia style discretization of gamma
epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, y)
  pmf[1] <- pgamma(1.5, alpha, rate = beta)
  for (i in 2:y) {
    pmf[i] <- pgamma(i+.5, alpha, rate = beta) - pgamma(i-.5, alpha, rate = beta)
  }
  
  pmf
}

zero_epidemia_gamma <- function(y, alpha, beta) {
  pmf <- rep(0, (y+1))
  pmf[1] <- pgamma(0.5, alpha, rate = beta)
  for (i in 2:(y+1)) {
    pmf[i] <- pgamma(i -1 +.5, alpha, rate = beta) - pgamma(i -1 -.5, alpha, rate = beta)
  }
  
  pmf
}


epidemia_hypoexp <- function(y, rates) {
  pmf <- rep(0, y)
  pmf[1] <- phypoexp(1.5, rates)
  for (i in 2:y) {
    pmf[i] <- phypoexp(i+.5, rates) - phypoexp(i-.5, rates)
  }
  
  pmf
}


epidemia_lognormal <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- plnorm(1.5, meanlog = params[1], sdlog = params[2])
  for (i in 2:y) {
    pmf[i] <- plnorm(i+.5, meanlog = params[1], sdlog = params[2]) - 
      plnorm(i-.5, meanlog = params[1], sdlog = params[2])
  }
  
  pmf
}

epidemia_weibull <- function(y, params) {
  pmf <- rep(0, y)
  pmf[1] <- pweibull(1.5, shape = params[1], scale = params[2])
  for (i in 2:y) {
    pmf[i] <- pweibull(i+.5, shape = params[1], scale = params[2]) - 
      pweibull(i-.5, shape = params[1], scale = params[2])
  }
  
  pmf
}



zero_epidemia_hypoexp <- function(y, rates) {
  pmf <- rep(0, (y+1))
  pmf[1] <- phypoexp(0.5, rates)
  for (i in 2:y+1) {
    pmf[i] <- phypoexp(i -1 +.5, rates) - phypoexp(i -1 -.5, rates)
  }
  
  pmf
}

