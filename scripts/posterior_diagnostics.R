library(tidyverse)
library(posterior)
library(kableExtra)

variable_name_key <-
  c(S_SEI = "$S_0$",
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
    α_t = "$\\exp\\left(\\tilde{\\alpha}_t\\right)$")

main_model_results <- 
  read_csv("results/generated_quantities/generated_quantities_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42_model_design=46_model_id=181_use_seroprev=true_use_tests=true.csv") %>% 
  rename_with(~str_c(".", .), c(iteration, chain)) %>% 
  as_draws()
 
summarized_results <-
  main_model_results %>% 
  summarize_draws("$\\hat{R}$" = rhat, "ESS" = ess_basic, .cores = 8) %>% 
  rename(Parameter = variable)
  
summarized_results %>%
  filter(str_detect(Parameter, "\\[\\d+\\]", negate = T)) %>% 
  mutate(Parameter = variable_name_key[Parameter]) %>% 
  kable(format = "latex", escape = F, caption = "example captions", booktabs = T, digits = 2) %>% 
  save_kable(file = "~/Documents/semi_parametric_COVID_19_OC_manuscript/tables/univariate_diagnostics.tex")


summarized_results %>% 
  rename(name_raw = Parameter) %>% 
  filter(str_detect(name_raw, "\\[\\d+\\]")) %>%
  mutate(parameter = name_raw %>% str_extract("^.+(?=\\[\\d+\\])"),
       index = name_raw%>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  select(-name_raw) %>% 
  filter(parameter %in% c("α_t", "Rₜ_t", "IFR_t")) %>% 
  mutate(parameter = variable_name_key[parameter]) %>% 
  pivot_longer(cols = c("$\\hat{R}$", "ESS")) %>%
  select(-index) %>% 
  group_by(parameter, name) %>% 
  summarize(min = min(value),
            avg =  mean(value),
            max = max(value)) %>%
  rename_with(str_to_title) %>% 
  rename_with(~str_c(., "."), .cols = c(Min, Avg, Max)) %>%
  pivot_wider(names_from = "Name", values_from = str_c(str_to_title(c("min", "avg", "max")), "."), names_sep = " ") %>% 
  relocate(Parameter, ends_with("$")) %>% 
  kable(format = "latex", escape = F, caption = "example captions", booktabs = T, digits = 2) %>% 
  save_kable(file = "~/Documents/semi_parametric_COVID_19_OC_manuscript/tables/multivariate_diagnostics.tex")
