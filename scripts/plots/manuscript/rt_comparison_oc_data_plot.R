library(tidyverse)
source("src/plot_functions.R")
model_table <- read_csv("model_table.csv")
max_t <- 42

dat <-
  read_csv("data/oc_data.csv") %>%
  select(time, cases, deaths, tests, date = end_date) %>%
  mutate(test_positivity = cases / tests)

dat_tidy <- 
  dat %>% 
  select(-date, -tests, -test_positivity) %>% 
  pivot_longer(-time)

time_interval_in_days <- as.numeric(dat$date[2] - dat$date[1])

time_date_key <- 
  dat %>% 
  select(time, date) %>% 
  add_row(., time = 0, date = min(.$date) - time_interval_in_days, .before = 1) 


posterior_generated_quantities_path <-
  tibble(full_path = dir_ls("results/tidy_posterior_generated_quantities")) %>%
  mutate(model_id = full_path %>%
    str_extract("(?<=model_id=)\\d+") %>%
    as.numeric()) %>%
  left_join(model_table %>% distinct(model_design, .keep_all = T)) %>%
  select(-model_id, -seed) %>%
  filter(
    constant_alpha == F,
    constant_IFR == F,
    constant_R0 == F,
    double_IFR_0 == F,
    half_alpha_0 == F,
    half_R0_0 == F,
    half_S_0 == F,
    max_t == 42,
    use_seroprev == T,
    use_tests == T
  ) %>%
  pull(full_path)

rt_intervals <-
  bind_rows(
    read_csv("results/rt_estim/rt_comparison_model_id=epiestim.csv") %>% 
      mutate(time = time - 1,
             name = "Rt",
             method = "EpiEstim",
             .width = 0.95) %>% 
      select(time, name, value = rt_median, .lower = rt_CI95l, .upper = rt_CI95u, .width, method) %>% 
      left_join(time_date_key),
    read_csv("results/rt_estim/rt_comparison_model_id=estimgamma.csv") %>% 
      mutate(time = time - 1) %>% 
      left_join(time_date_key),
    read_csv("results/rt_estim/rt_comparison_model_id=estimnormal.csv") %>% 
      mutate(time = time - 1) %>% 
      left_join(time_date_key),
    read_csv(posterior_generated_quantities_path) %>%
      filter(name == "Rₜ_t") %>%
      select(name, date, value, .lower, .upper, .width) %>%
      mutate(method = "Full Model") %>% 
      left_join(time_date_key) %>% 
      filter(time <= 41)
  ) %>% 
  mutate(method = method |> 
           fct_inorder() |> 
           fct_recode(`Rt-estim-gamma` = "estim_gamma",
                      Epidemia = "estim_normal"))

rt_comparison_oc_data_plot <- 
  rt_intervals %>% 
  filter(method != "Rt-estim-gamma", time > 0) %>% 
  filter(.width == 0.95) %>% 
  ggplot(aes(x = date, y = value, ymin = .lower, ymax = .upper, fill = method)) +
  geom_lineribbon(alpha = 0.5) +
  scale_y_continuous(name = my_labeller["Rₜ_t"]) +
  ggtitle("Posterior Effective Reproduction Number",
          subtitle = "95% Credible Intervals") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_date(name = "Date") +
  scale_fill_discrete(name = "Method", label = str_to_title) +
  theme(legend.position = "bottom")

save_plot(filename = path(figures_dir, "rt_comparison_oc_data_plot", ext = "pdf"),
          plot = rt_comparison_oc_data_plot)
