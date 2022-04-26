library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)

options(mc.cores = parallel::detectCores())
dat <- read_csv("data/oc_data.csv") %>% 
  filter(time <= 42)

desired_quantiles <- c(0.05, 0.5, 0.95)



# Deaths ------------------------------------------------------------------
deaths_model_nb <- brm(bf(deaths ~ s(time)),
                       data = dat,
                       family = negbinomial(),
                       cores = 4,
                       seed = 17,
                       iter = 4000,
                       warmup = 1000,
                       control = list(adapt_delta = 0.99, max_treedepth = 12),
                       backend = "cmdstanr")


deaths_model_nb_shape_draws <- gather_draws(deaths_model_nb, shape)$.value

deaths_model_nb_mean <- mean(log(deaths_model_nb_shape_draws))
deaths_model_nb_sd <- sd(log(deaths_model_nb_shape_draws))

tibble(desired = quantile(deaths_model_nb_shape_draws, desired_quantiles),
       obtained = exp(qnorm(p = desired_quantiles, mean = deaths_model_nb_mean, deaths_model_nb_sd)))



# Cases Beta-Binomial -----------------------------------------------------
# https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html

beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)

stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"

posterior_predict_beta_binomial2 <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  trials <- prep$data$vint1[i]
  beta_binomial2_rng(mu, phi, trials)
}

stanvars <- stanvar(scode = stan_funs, block = "functions")

cases_model_bb <- brm(
  cases | vint(tests) ~ s(time), data = dat, 
  family = beta_binomial2, stanvars = stanvars,
  cores = 4,
  seed = 18,
  iter = 8000,
  warmup = 4000,
  control = list(adapt_delta = 0.8)
)

expose_functions(cases_model_bb, vectorize = TRUE)

log_lik_beta_binomial2 <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, trials)
}

cases_model_bb_phi_draws <- gather_draws(cases_model_bb, phi)$.value

cases_model_bb_mean <- mean(log(cases_model_bb_phi_draws))
dput(cases_model_bb_mean)

cases_model_bb_sd <- sd(log(cases_model_bb_phi_draws))
dput(cases_model_bb_sd)

tibble(desired = quantile(cases_model_bb_phi_draws, desired_quantiles),
       obtained = exp(qnorm(p = desired_quantiles, mean = cases_model_bb_mean, cases_model_bb_sd)))

# Cases - Negative-Binomial -----------------------------------------------
cases_model_nb <- brm(bf(cases ~ s(time)),
                      data = dat,
                      family = negbinomial(),
                      cores = 4,
                      seed = 17,
                      iter = 4000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.99, max_treedepth = 12),
                      backend = "cmdstanr")





cases_model_nb_shape_draws <- gather_draws(cases_model_nb, shape)$.value

cases_model_nb_mean <- mean(log(cases_model_nb_shape_draws))

cases_model_nb_sd <- sd(log(cases_model_nb_shape_draws))


tibble(desired = quantile(cases_model_nb_shape_draws, desired_quantiles),
       obtained = exp(qnorm(p = desired_quantiles, mean = cases_model_nb_mean, cases_model_nb_sd)))



# Plot Posterior Predictive -----------------------------------------------
predicted_draws <- 
  predicted_draws(deaths_model_nb, dat, value = "deaths_nb") %>% 
  left_join(predicted_draws(cases_model_bb, dat, value = "cases_bb")) %>% 
  left_join(predicted_draws(cases_model_nb, dat, value = "cases_nb")) %>% 
  ungroup() %>% 
  mutate(across(starts_with("cases"), ~{. / tests})) %>% 
  rename_all(~str_replace(., "cases", "pos")) %>% 
  select(date = end_date, starts_with("pos"), starts_with("deaths"))


predicted_draws_data <- 
  predicted_draws %>% 
  select(-contains("_")) %>% 
  distinct() %>% 
  mutate(pos_nb = pos,
         pos_bb = pos,
         deaths_nb = deaths) %>% 
  select(-pos, -deaths) %>% 
  pivot_longer(-date,
               names_to = c("data", "model"),
               names_sep = "_")

predicted_draws_intervals <- 
  predicted_draws %>% 
  select(date, contains("_")) %>% 
  pivot_longer(-date,
               names_to = c("data", "model"),
               names_sep = "_") %>% 
  group_by(date, data, model) %>% 
  median_qi(.width = c(0.5, 0.8, 0.95))

ggplot(mapping = aes(date, value)) +
  facet_wrap(data ~ model, scales = "free_y") +
  geom_lineribbon(data = predicted_draws_intervals, mapping = aes(ymin = .lower, ymax = .upper)) +
  geom_point(data = predicted_draws_data) +
  scale_fill_brewer() +
  cowplot::theme_cowplot()

# Save Results ------------------------------------------------------------

# maybe write a .jl file?
tibble(model = c("cases_bb", "cases_nb", "deaths_nb"),
       mean = c(cases_model_nb_mean, deaths_model_nb_mean, cases_model_bb_mean),
       sd = c(cases_model_nb_sd, deaths_model_nb_sd, cases_model_bb_sd))

c(str_c("const ϕ_cases_nb_non_centered_mean = ", cases_model_nb_mean),
  str_c("const ϕ_cases_nb_non_centered_sd = ", cases_model_nb_sd),
  str_c("const ϕ_cases_bb_non_centered_mean = ", cases_model_bb_mean),
  str_c("const ϕ_cases_bb_non_centered_sd = ", cases_model_bb_sd),
  str_c("const ϕ_deaths_non_centered_mean = ", deaths_model_nb_mean),
  str_c("const ϕ_deaths_non_centered_sd = ", deaths_model_nb_sd)) %>% 
  write_lines("src/prior_constants_overdispersion.jl")
