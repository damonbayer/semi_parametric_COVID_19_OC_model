library(tidyverse)
source("src/plot_functions.R")
library(posterior)
library(bayesplot)
library(tidybayes)
library(gridExtra)

file_path_posterior_gq <- "experiments/fixed_sigma_experiment/results/posterior_generated_quantities/posterior_generated_quantities_fixed_sigma_all_tight_init.csv"
file_path_prior_gq <- "experiments/fixed_sigma_experiment/results/prior_generated_quantities/prior_generated_quantities_fixed_sigma_all_tight_init.csv"
file_path_posterior_predictive <- "experiments/fixed_sigma_experiment/results/posterior_predictive/posterior_predictive_fixed_sigma_all_tight_init.csv"
figure_path <- path("experiments/fixed_sigma_experiment/results/comprehensive_model_diagnostics_fixed_sigma_all_tight_init", ext = "pdf")

dat <- 
  read_csv("data/oc_data.csv") %>% 
  filter(time <= 42) %>% 
  select(time, cases, deaths, tests) %>% 
  mutate(test_positivity = cases / tests)

# for simulated data
# dat <- 
#   read_csv("experiments/fixed_sigma_experiment/simulated_data.csv") %>% 
#   select(-iteration, -chain) %>% 
#   pivot_longer(everything()) %>% 
#   separate(col = name,
#            into = c("name", "index"),
#            sep = "\\[|\\]",
#            remove = T,
#            fill = "right",
#            extra = "drop",
#            convert = T) %>% 
#   mutate(name = str_remove(name, "data_new_")) %>% 
#   rename(time = index) %>% 
#   pivot_wider(names_from = name, values_from = value) %>% 
#   select(-data_seroprev_cases) %>% 
#   left_join(read_csv("data/oc_data.csv") %>% select(time, tests)) %>% 
#   mutate(test_positivity = cases / tests)

fix_gq_names <- function(x) {
  x %>% 
    str_replace("₀", "0") %>% 
    str_replace("ₜ", "t") %>% 
    str_replace("α", "alpha") %>% 
    str_replace("β", "beta") %>% 
    str_replace("ρ", "rho") %>% 
    str_replace("σ", "sigma") %>% 
    str_replace("ϕ", "phi")
}

create_gq_draws_for_plotting <- function(file_path) {
  read_csv(file_path) |> 
    rename_with(~str_c(".", .), c(iteration, chain)) |> 
    as_draws() %>% 
    rename_with(fix_gq_names)
}

create_predictive_draws_for_plotting <- function(file_path) {
  read_csv(file_path) |> 
    rename_with(~str_c(".", .), c(iteration, chain)) |> 
    as_draws()
}

tidy_gq_by_chain <- function(gq_draws_for_plotting) {
  gq_draws_for_plotting |> 
    pivot_longer(-starts_with(".")) |> 
    group_by(name, .chain) |> 
    median_qi(.width = 0.8) |>
    separate(col = name,
             into = c("name", "index"),
             sep = "\\[|\\]",
             remove = T,
             fill = "right",
             extra = "drop",
             convert = T) |> 
    mutate(.chain = factor(.chain)) |> 
    arrange(.chain, index)
}

tidy_predictive_by_chain <- function(predictive_draws_for_plotting) {
  predictive_draws_for_plotting %>% 
    select(-matches("data_seroprev_cases\\[\\d+\\]")) %>% 
    pivot_longer(-starts_with(".")) %>% 
    separate(col = name,
             into = c("name", "index"),
             sep = "\\[|\\]",
             remove = T,
             fill = "right",
             extra = "drop",
             convert = T) %>% 
    mutate(name = str_remove(name, "data_new_")) %>% 
    rename(time = index) %>% 
    pivot_wider(names_from = name, values_from = value) %>% 
    select(-c(.draw, .iteration)) %>% 
    left_join(select(dat, time, tests)) %>% 
    mutate(test_positivity = cases / tests) %>% 
    select(-cases, -tests) %>% 
    pivot_longer(cols = c(deaths, test_positivity)) %>% 
    group_by(.chain, time, name) %>% 
    median_qi(.width = 0.8) %>% 
    mutate(.chain = factor(.chain))
}

dat_tidy <- 
  dat %>% 
  select(-cases, -tests) %>% 
  pivot_longer(-time)

gq_posterior_draws_for_plotting <- create_gq_draws_for_plotting(file_path_posterior_gq)
gq_prior_draws_for_plotting <- create_gq_draws_for_plotting(file_path_prior_gq)

gq_posterior_tidy_by_chain <- tidy_gq_by_chain(gq_posterior_draws_for_plotting)
gq_prior_tidy_by_chain <- tidy_gq_by_chain(gq_prior_draws_for_plotting)

posterior_predictive_draws_for_plotting <- create_predictive_draws_for_plotting(file_path_posterior_predictive)
posterior_predictvie_tidy_by_chain <- tidy_predictive_by_chain(posterior_predictive_draws_for_plotting)


# Diagnostic Summaries ----------------------------------------------------
summarized_results <-
  gq_posterior_draws_for_plotting %>%
  summarize_draws(rhat, ess_basic, .cores = 8)

univariate_ess_rhat <- 
  summarized_results |> 
  filter(str_detect(variable, "\\[\\d+\\]", negate = T)) %>% 
  drop_na()

multivariate_ess_rhat <- 
  summarized_results %>%
  rename(name_raw = variable) %>%
  filter(str_detect(name_raw, "\\[\\d+\\]")) %>%
  mutate(
    parameter = name_raw %>% str_extract("^.+(?=\\[\\d+\\])"),
    index = name_raw %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric()
  ) %>%
  select(-name_raw) %>%
  filter(parameter %in% c("alpha_t", "Rt_t", "IFR_t")) |> 
  pivot_longer(cols = c(rhat, ess_basic)) %>%
  select(-index) %>%
  group_by(parameter, name) %>%
  summarize(
    min = min(value),
    avg = mean(value),
    med = median(value),
    max = max(value)
  ) |> 
  arrange(name, parameter) %>% 
  ungroup()

# Plots -------------------------------------------------------------------
univariate_gq_names <-
  tibble(names = colnames(gq_posterior_draws_for_plotting)) %>% 
  filter(str_starts(names, "\\.", negate = T)) %>% 
  filter(str_detect(names, "[^\\_t]\\[\\d+\\]", negate = T)) %>% 
  filter(str_detect(names, "\\[\\d+\\]", negate = T) | str_detect(names, "\\[0\\]")) %>% 
  distinct(names) %>% 
  pull(names)

univariate_gq_tidy_names <- 
  tibble(names = colnames(gq_posterior_draws_for_plotting)) %>% 
  filter(str_starts(names, "\\.", negate = T)) %>% 
  filter(str_detect(names, "[^\\_t]\\[\\d+\\]", negate = T)) %>% 
  mutate(names = str_remove_all(names, "\\[\\d+\\]")) %>% 
  distinct(names) %>% 
  pull(names)

compartment_names <- c("S", "E", "I", "R", "D")

univariate_gq_trace_plot <- 
  gq_posterior_draws_for_plotting |> 
  select(starts_with("."), all_of(univariate_gq_names)) |> 
  thin_draws(5) |>
  tidy_draws() |> 
  pivot_longer(-starts_with(".")) |> 
  mutate(.chain = factor(.chain)) |> 
  ggplot(aes(.iteration, value, color = .chain)) +
  facet_grid(name ~ .chain, scales = "free_y") +
  geom_line() +
  theme_cowplot() +
  ggtitle("Univariate GQ Trace Plot")

univariate_interval_plot <- 
  bind_rows(gq_posterior_tidy_by_chain,
            gq_prior_tidy_by_chain %>% 
              mutate(.chain = "prior")) %>% 
  filter(name %in% univariate_gq_tidy_names) %>% 
  mutate(name = if_else(str_ends(name, "_t"), str_c(name, "_init"), name)) %>% 
  ggplot(aes(value, .chain, xmin = .lower, xmax =.upper, color = .chain)) +
  facet_wrap(~name, scales = "free_x") +
  geom_interval() +
  ggtitle("Univariate Eye Plot")

univariate_pairs_plot <- 
  gq_posterior_draws_for_plotting |> 
  select(starts_with("."), all_of(univariate_gq_names)) %>% 
  mcmc_pairs(off_diag_fun = "hex")

compartments_by_chain_plot <- 
  gq_posterior_tidy_by_chain |> 
  filter(!is.na(index)) |> 
  filter(name %in% compartment_names) |> 
  mutate(name = name |> fct_relevel(compartment_names)) |> 
  ggplot(aes(index, value, ymin = .lower, ymax = .upper, fill = .chain)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon(alpha = 0.5) +
  scale_y_continuous(labels = comma) +
  ggtitle("Compartments by Chain")

time_varying_parameters_by_chain_plot <- 
  gq_posterior_tidy_by_chain |> 
  filter(!is.na(index)) |> 
  filter(str_ends(name, "_t")) |> 
  ggplot(aes(index, value, ymin = .lower, ymax = .upper, fill = .chain)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon(alpha = 0.5) +
  scale_y_continuous(labels = comma) +
  ggtitle("Time-Varying Parameters by Chain")

univariate_acf_plot <-
  gq_posterior_draws_for_plotting |> 
  select(starts_with("."), all_of(univariate_gq_names)) %>% 
  mcmc_acf(x = ., lags = max(.$.iteration) - 2) + 
  ggtitle("Univariate ACF")

univariate_ess_rhat_plot <- 
  univariate_ess_rhat %>% 
  pivot_longer(-variable, names_to = "value_type") %>% 
  ggplot(aes(variable, value)) +
  facet_wrap(~value_type, scales = "free_y", ncol = 1)+
  geom_col() +
  ggtitle("Univariate ESS rhat")

multivariate_ess_rhat_plot <- 
  multivariate_ess_rhat %>% 
  pivot_longer(-c(parameter, name), names_to = "value_type") %>% 
  ggplot(aes(value, parameter, color = value_type)) +
  facet_wrap(~name, scales = "free_x", ncol = 1) +
  geom_point(size = 5) +
  ggtitle("Multivariate ESS rhat")

posterior_predictive_plot <- 
  ggplot(mapping = aes(time, value)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon(data = posterior_predictvie_tidy_by_chain,
                  mapping = aes(ymin = .lower, ymax = .upper, fill = .chain), alpha = 0.25, step = "mid") +
  geom_point(data = dat_tidy) +
  scale_y_continuous() +
  ggtitle("Posterior Predictive")


# Save Plots --------------------------------------------------------------
ggsave2(filename = figure_path,
        plot = ls()[str_ends(ls(), "_plot")] |> 
          map(get) |>
          marrangeGrob(ncol = 1, nrow = 1),
        width = 16*1.5,
        height = 9*1.5)
