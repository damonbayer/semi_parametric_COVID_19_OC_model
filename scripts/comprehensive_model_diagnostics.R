source("src/plot_functions.R")
library(tidyverse)
library(posterior)
library(bayesplot)
library(tidybayes)
library(gridExtra)

file_path <- "experiments/blocking_strategies_experiment/results/gq/gq_init_ess_sigma_hmc.csv"

draws_for_plotting <- 
  read_csv(file_path) |> 
  rename_with(~str_c(".", .), c(iteration, chain)) |> 
  as_draws()

draws_for_plotting_tidy_by_chain <-
  draws_for_plotting |> 
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

univariate_trace_plot <- 
  draws_for_plotting |> 
  select(-matches("\\[\\d+\\]")) |> 
  thin_draws(5) |>
  tidy_draws() |> 
  pivot_longer(-starts_with(".")) |> 
  mutate(.chain = factor(.chain)) |> 
  ggplot(aes(.iteration, value, color = .chain)) +
  facet_grid(name ~ .chain, scales = "free_y") +
  geom_line() +
  theme_cowplot()

univariate_eye_plot <- 
  draws_for_plotting |> 
  select(-matches("\\[\\d+\\]"), `R₀_t[1]`) |> 
  tidy_draws() |> 
  pivot_longer(-starts_with(".")) |> 
  mutate(.chain = factor(.chain)) |> 
  ggplot(aes(value, fill = .chain, color = .chain)) +
  facet_wrap(~ name, scales = "free_x") +
  stat_halfeye(alpha = 0.5, normalize = "panels") +
  theme_cowplot()

univariate_pairs_plot <- 
  draws_for_plotting |> 
  select(-matches("\\[\\d+\\]"), `R₀_t[1]`) |> 
  mcmc_pairs(off_diag_fun = "hex")

compartments_by_chain_plot <- 
  draws_for_plotting_tidy_by_chain |> 
  filter(!is.na(index)) |> 
  filter(name %in% c("S", "E", "I", "R", "D")) |> 
  mutate(name = name |> fct_relevel(c("S", "E", "I", "R", "D"))) |> 
  ggplot(aes(index, value, ymin = .lower, ymax = .upper, fill = .chain)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon(alpha = 0.5) +
  scale_y_continuous(labels = comma)

time_varying_parameters_by_chain_plot <- 
  draws_for_plotting_tidy_by_chain |> 
  filter(!is.na(index)) |> 
  filter(str_ends(name, "_t")) |> 
  ggplot(aes(index, value, ymin = .lower, ymax = .upper, fill = .chain)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon(alpha = 0.5) +
  scale_y_continuous(labels = comma)

univariate_acf_plot <-
  draws_for_plotting |> 
  thin_draws(10) |> 
  select(-matches("\\[\\d+\\]")) %>% 
  mcmc_acf(x = ., lags = max(.$.iteration) - 2)

summarized_results <-
  draws_for_plotting %>%
  summarize_draws(rhat, ess_basic, .cores = 8)

univariate_ess_rhat <- 
  summarized_results |> 
  filter(str_detect(variable, "\\[\\d+\\]", negate = T))

multivariate__ess_rhat <- 
  summarized_results %>%
  rename(name_raw = variable) %>%
  filter(str_detect(name_raw, "\\[\\d+\\]")) %>%
  mutate(
    parameter = name_raw %>% str_extract("^.+(?=\\[\\d+\\])"),
    index = name_raw %>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric()
  ) %>%
  select(-name_raw) %>%
  filter(parameter %in% c("α_t", "Rₜ_t", "IFR_t")) |> 
  pivot_longer(cols = c(rhat, ess_basic)) %>%
  select(-index) %>%
  group_by(parameter, name) %>%
  summarize(
    min = min(value),
    avg = mean(value),
    med = median(value),
    max = max(value)
  ) |> 
  arrange(name, parameter)


ggsave2(filename = path(str_c(path_ext_remove(file_path), "_comprehensive_model_diagnostics"), ext = "pdf"),
        plot = ls()[str_ends(ls(), "_plot")] |> 
          map(get) |>
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 12)
