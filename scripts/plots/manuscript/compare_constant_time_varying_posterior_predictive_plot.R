library(tidyverse)
library(tidybayes)
library(scales)
source("src/plot_functions.R")

max_t <- 42

dat <- 
  read_csv("data/oc_data.csv") %>% 
  filter(time <= max_t) %>% 
  mutate(index = 1:n()) %>% 
  rename(date = end_date) %>% 
  select(-start_date)

max_date <- max(dat$date)

model_table <-
  read_csv("model_table.csv") %>% 
  select(-model_id, -seed) %>% 
  distinct()

dat_tidy <-
  dat %>% 
  select(date, cases, tests, deaths) %>% 
  mutate(test_positivity = cases / tests) %>% 
  select(-tests) %>% 
  pivot_longer(-date)

all_posterior_predictive <- 
  tibble(full_path = dir_ls("results/tidy_posterior_predictive")) %>% 
  mutate(model_design = full_path %>% 
           str_extract("(?<=model_design=)\\d+") %>% 
           as.integer()) %>% 
  left_join(model_table) %>% 
  filter(max_t == 42,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F,
         use_seroprev == T,
         use_tests == T) %>% 
  select(full_path, model_design, starts_with("constant")) %>% 
  mutate("all_time_varying" = !constant_alpha & !constant_IFR & !constant_R0) %>% 
  pivot_longer(-c(full_path, model_design), names_to = "model_name") %>% 
  filter(value) %>% 
  select(-value) %>% 
  mutate(posterior_predictive = map(full_path, read_csv)) %>% 
  select(-full_path) %>% 
  unnest(posterior_predictive) %>% 
  filter(date <= max_date) %>% 
  bind_rows(
    .,
    filter(., name == "cases") %>%
      mutate(name = "test_positivity") %>%
      left_join(dat %>% select(date, tests)) %>%
      mutate(across(c(value, .lower, .upper), ~ . / tests)) %>%
      select(-tests)
  )

deaths_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = all_posterior_predictive %>% 
                    filter(name == "deaths",
                           .width == 0.8),
                  mapping = aes(ymin = .lower, ymax = .upper,
                                fill = model_name, color = model_name),
                  alpha = 0.25, key_glyph = "rect") +
  geom_point(data = dat_tidy %>% filter(name == "deaths")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date") +
  scale_fill_discrete(name = "Model",
                      labels = TeX(c(all_time_varying = "All Time-Varying",
                                     constant_R0 = "Constant $R_0$",
                                     constant_alpha = "Constant $\\alpha",
                                     constant_IFR  = "Constant $\\eta$"))) +
  scale_color_discrete(name = "Model",
                       labels = TeX(c(all_time_varying = "All Time-Varying",
                                      constant_R0 = "Constant $R_0$",
                                      constant_alpha = "Constant $\\alpha",
                                      constant_IFR  = "Constant $\\eta$"))) +
  ggtitle("Deaths - Sensitivity Analysis", subtitle = "80% Posterior Predictive Credible Intervals")


test_positivity_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = all_posterior_predictive %>% 
                    filter(name == "test_positivity",
                           .width == 0.8),
                  mapping = aes(ymin = .lower, ymax = .upper,
                                fill = model_name, color = model_name),
                  alpha = 0.25, key_glyph = "rect") +
  geom_point(data = dat_tidy %>% filter(name == "test_positivity")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(name = "Test Positivity", labels = percent) +
  scale_x_date(name = "Date") +
  scale_fill_discrete(name = "Model",
                      labels = TeX(c(all_time_varying = "All Time-Varying",
                                     constant_R0 = "Constant $R_0$",
                                     constant_alpha = "Constant $\\alpha",
                                     constant_IFR  = "Constant $\\eta$"))) +
  scale_color_discrete(name = "Model",
                       labels = TeX(c(all_time_varying = "All Time-Varying",
                                      constant_R0 = "Constant $R_0$",
                                      constant_alpha = "Constant $\\alpha",
                                      constant_IFR  = "Constant $\\eta$"))) +
  ggtitle("Test Positivity - Sensitivity Analysis", subtitle = "80% Posterior Predictive Credible Intervals") +
  theme(legend.position = "none")


compare_constant_time_varying_posterior_predictive_plot <- 
  plot_grid(test_positivity_plot, deaths_plot, align = "hv", nrow = 2, ncol = 1)

save_plot(filename = path(figures_dir, "compare_constant_time_varying_posterior_predictive_plot", ext = "pdf"),
          plot = compare_constant_time_varying_posterior_predictive_plot,
          ncol = 1, nrow = 2, base_asp = 2)
