library(tidyverse)
source("src/plot_functions.R")
library(posterior)
library(bayesplot)
library(tidybayes)
library(gridExtra)

file_path_gq <- "experiments/fixed_sigma_experiment/results/gq/gq_fixed_sigma_all_tight_init.csv"
file_path_pp <- "experiments/fixed_sigma_experiment/results/pp/pp_fixed_sigma_all_tight_init.csv"


dat <- 
  read_csv("data/oc_data.csv") %>% 
  filter(time <= 42) %>% 
  select(time, cases, deaths, tests) %>% 
  mutate(test_positivity = cases / tests)

dat_tidy <- 
  dat %>% 
  select(-cases, -tests) %>% 
  pivot_longer(-time)

gq_draws_for_plotting <- 
  read_csv(file_path_gq) |> 
  rename_with(~str_c(".", .), c(iteration, chain)) |> 
  as_draws() %>% 
  rename_with(~ str_replace(., "₀", "0") %>% 
                str_replace("ₜ", "t") %>% 
                str_replace("α", "alpha") %>% 
                str_replace("β", "beta") %>% 
                str_replace("ρ", "rho") %>% 
                str_replace("σ", "sigma") %>% 
                str_replace("ϕ", "phi"))

gq_draws_for_plotting_tidy_by_chain <-
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

pp_draws_for_plotting <- 
  read_csv(file_path_pp) |> 
  rename_with(~str_c(".", .), c(iteration, chain)) |> 
  as_draws()


pp_draws_for_plotting_tidy_by_chain <- 
  pp_draws_for_plotting %>% 
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


univariate_trace_plot <- 
  gq_draws_for_plotting |> 
  select(-matches("\\[\\d+\\]"), `R0_t[0]`, `IFR_t[0]`, `alpha_t[0]`) |> 
  thin_draws(5) |>
  tidy_draws() |> 
  pivot_longer(-starts_with(".")) |> 
  mutate(.chain = factor(.chain)) |> 
  ggplot(aes(.iteration, value, color = .chain)) +
  facet_grid(name ~ .chain, scales = "free_y") +
  geom_line() +
  theme_cowplot() +
  ggtitle("Univariate Trace Plot")

univariate_eye_plot <- 
  gq_draws_for_plotting |> 
  select(-matches("\\[\\d+\\]"), `R0_t[0]`, `IFR_t[0]`, `alpha_t[0]`) |> 
  tidy_draws() |> 
  pivot_longer(-starts_with(".")) |> 
  mutate(.chain = factor(.chain)) |> 
  ggplot(aes(value, fill = .chain, color = .chain)) +
  facet_wrap(~ name, scales = "free_x") +
  stat_halfeye(alpha = 0.5, normalize = "panels") +
  theme_cowplot() +
  ggtitle("Univariate Eye Plot")

univariate_pairs_plot <- 
  gq_draws_for_plotting |> 
  select(-matches("\\[\\d+\\]"), `R0_t[1]`) |> 
  mcmc_pairs(off_diag_fun = "hex")

compartments_by_chain_plot <- 
  gq_draws_for_plotting_tidy_by_chain |> 
  filter(!is.na(index)) |> 
  filter(name %in% c("S", "E", "I", "R", "D")) |> 
  mutate(name = name |> fct_relevel(c("S", "E", "I", "R", "D"))) |> 
  ggplot(aes(index, value, ymin = .lower, ymax = .upper, fill = .chain)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon(alpha = 0.5) +
  scale_y_continuous(labels = comma) +
  ggtitle("Compartments by Chain")

time_varying_parameters_by_chain_plot <- 
  gq_draws_for_plotting_tidy_by_chain |> 
  filter(!is.na(index)) |> 
  filter(str_ends(name, "_t")) |> 
  ggplot(aes(index, value, ymin = .lower, ymax = .upper, fill = .chain)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon(alpha = 0.5) +
  scale_y_continuous(labels = comma) +
  ggtitle("Time-Varying Parameters by Chain")

univariate_acf_plot <-
  gq_draws_for_plotting |> 
  select(-matches("\\[\\d+\\]")) %>% 
  mcmc_acf(x = ., lags = max(.$.iteration) - 2)

summarized_results <-
  gq_draws_for_plotting %>%
  summarize_draws(rhat, ess_basic, .cores = 8)

summarized_results |> 
  filter(variable == "R0_t[1]")

univariate_ess_rhat <- 
  summarized_results |> 
  filter(str_detect(variable, "\\[\\d+\\]", negate = T)) %>% 
  drop_na()

univariate_ess_rhat_plot <- 
  univariate_ess_rhat %>% 
  pivot_longer(-variable, names_to = "value_type") %>% 
  ggplot(aes(variable, value)) +
  facet_wrap(~value_type, scales = "free_y", ncol = 1)+
  geom_col() +
  ggtitle("Univariate ESS rhat")

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

multivariate_ess_rhat_plot <- 
  multivariate_ess_rhat %>% 
  pivot_longer(-c(parameter, name), names_to = "value_type") %>% 
  ggplot(aes(value, parameter, color = value_type)) +
  facet_wrap(~name, scales = "free_x", ncol = 1) +
  geom_point(size = 5) +
  ggtitle("Multivariate ESS rhat")

pp_plot <- 
  ggplot(mapping = aes(time, value)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon(data = pp_draws_for_plotting_tidy_by_chain,
                  mapping = aes(ymin = .lower, ymax = .upper, fill = .chain), alpha = 0.25, step = "mid") +
  geom_point(data = dat_tidy) +
  scale_y_continuous() +
  ggtitle("Posterior Predictive")

ls()[str_ends(ls(), "_plot")]

ggsave2(filename = path(str_c(path_ext_remove(file_path_gq), "_comprehensive_model_diagnostics"), ext = "pdf"),
        plot = ls()[str_ends(ls(), "_plot")] |> 
          map(get) |>
          marrangeGrob(ncol = 1, nrow = 1),
        width = 16,
        height = 9)
