library(tidyverse)
library(posterior)
library(kableExtra)
source("src/plot_functions.R")

target_model_design <- 46

variable_name_key <-
  c(
    S_SEI = "$S_0$",
    I_EI = "$\\tilde{I}_{0}$",
    dur_latent_days = "$1 / \\gamma$",
    dur_infectious_days = "$1 / \\nu$",
    ϕ_deaths = "$\\phi_D$",
    ρ_death = "$\\rho^D$",
    σ_α = "$\\sigma_\\alpha$",
    ϕ_cases_bb = "$\\phi_C$",
    σ_R0 = "$\\sigma_{R_0}$",
    σ_IFR = "$\\sigma_\\eta$",
    Rₜ_t = "$\\exp\\left(\\tilde{R}_{t,t}\\right)$",
    IFR_t = "$\\expit\\left(\\tilde{\\eta}_t\\right)$",
    α_t = "$\\exp\\left(\\tilde{\\alpha}_t\\right)$"
  )

posterior_lp <- 
  dir_ls("results/posterior_lp") %>% 
  enframe(name = NULL) %>% 
  filter(str_detect(value, str_c("model_design=", target_model_design, "_"))) %>% 
  pull(value) %>% 
  read_csv() %>% 
  as_draws()

posterior_diagnostics <- 
  enframe(dir_ls("results/posterior_diagnostics"), name = NULL) %>% 
  filter(str_detect(value, str_c("model_design=", target_model_design, "_"))) %>% 
  pull(value) %>% 
  read_csv()

posterior_diagnostics %>% 
  filter(is.na(date)) %>% 
  filter(name != "lp") %>% 
  mutate(Parameter = variable_name_key[name]) %>% 
  select(Parameter, "$\\hat{R}$" = rhat, "ESS" = ess_basic) %>% 
  kable(
    format = "latex",
    escape = F,
    caption = "Convergence diagnostics for time-stationary parameters for the main model fit to the Orange County data set.",
    booktabs = T,
    digits = 2,
    label = "univariate_diagnostics"
  ) %>%
  save_kable(file = "~/Documents/semi_parametric_COVID_19_OC_manuscript/tables/univariate_diagnostics.tex")

posterior_diagnostics %>% 
  filter(name %in% c("α_t", "Rₜ_t", "IFR_t")) %>% 
  mutate(Parameter = variable_name_key[name]) %>% 
  select(Parameter, "$\\hat{R}$" = rhat, "ESS" = ess_basic) %>% 
  drop_na() %>% 
  group_by(Parameter) %>% 
  summarize(across(where(is.numeric), list(min = min,  avg = mean, max = max), .names = "{str_to_title(.fn)}. {.col}")) %>% 
  kable(format = "latex",
        escape = F,
        caption = "Convergence diagnostics for time-varying parameters for the main model fit to the Orange County data set.",
        booktabs = T,
        digits = 2,
        label = "multivariate_diagnostics") %>%
  save_kable(file = "~/Documents/semi_parametric_COVID_19_OC_manuscript/tables/multivariate_diagnostics.tex")
