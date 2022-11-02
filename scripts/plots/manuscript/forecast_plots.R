library(tidyverse)
source("src/plot_functions.R")

model_table <- 
  read_csv("model_table.csv") %>% 
  select(-model_id, -seed) %>% 
  distinct()

time_date_converter <- 
  read_csv("data/oc_data.csv") %>% 
  select(time, date = end_date)

dat_tidy <-
  read_csv("data/oc_data.csv") %>% 
  rename(date = end_date) %>% 
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
  filter(constant_alpha == F,
         constant_IFR == F,
         constant_R0 == F,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F,
         use_seroprev == T,
         use_tests == T) %>% 
  mutate(pp = map2(full_path, max_t,
                   ~read_csv(.x) %>% 
                       left_join(time_date_converter) %>% 
                       mutate(weeks_ahead = time - .y) %>% 
                       filter(weeks_ahead > 0))) %>% 
  select(model_design, max_t, pp) %>% 
  unnest(pp)

weeks_ahead_dat <- 
  all_posterior_predictive %>% 
  group_by(weeks_ahead) %>% 
  summarize(min_date = min(date),
            max_date = max(date)) %>% 
  crossing(dat_tidy) %>% 
  filter(date >= min_date,
         date <= max_date) %>% 
  select(-ends_with("_date"))

weeks_ahead_to_plot <- c(1, 4)

deaths_forecast_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = all_posterior_predictive %>% 
                    filter(name == "deaths",
                           weeks_ahead %in% weeks_ahead_to_plot),
                  mapping = aes(ymin = .lower, ymax = .upper),
                  step = "mid") +
  geom_point(data = weeks_ahead_dat %>%
               filter(name == "deaths",
                      weeks_ahead %in% weeks_ahead_to_plot)) +
  facet_wrap(~ weeks_ahead, scales = "free_y",
             labeller = as_labeller(~str_c(., " Week", if_else(. == 1, "", "s"), " Ahead Forecasts")),
             ncol = 1) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date") +
  ggtitle("Deaths Forecasts") +
  my_theme


test_positivity_forecast_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = all_posterior_predictive %>% 
                    filter(name == "test_positivity",
                           weeks_ahead %in% weeks_ahead_to_plot),
                  mapping = aes(ymin = .lower, ymax = .upper),
                  step = "mid") +
  geom_point(data = weeks_ahead_dat %>%
               filter(name == "test_positivity",
                      weeks_ahead %in% weeks_ahead_to_plot)) +
  facet_wrap(~ weeks_ahead, scales = "free_y",
             labeller = as_labeller(~str_c(., " Week", if_else(. == 1, "", "s"), " Ahead Forecasts")),
             ncol = 1) +
  scale_y_continuous(name = "Test Positivity", labels = percent) +
  scale_x_date(name = "Date") +
  ggtitle("Test Positivity Forecasts") +
  my_theme

combined_forecast_plot <- plot_grid(deaths_forecast_plot,
                                    test_positivity_forecast_plot,
                                    ncol = 2,
                                    align = "hv")

save_plot(
  filename = path(figures_dir, "deaths_forecast_plot", ext = "pdf"),
  plot = deaths_forecast_plot,
  ncol = 1,
  nrow = 2)

save_plot(
  filename = path(figures_dir, "test_positivity_forecast_plot", ext = "pdf"),
  plot = test_positivity_forecast_plot,
  ncol = 1,
  nrow = 2)

save_plot(
  filename = path(figures_dir, "combined_forecast_plot", ext = "pdf"),
  plot = combined_forecast_plot,
  ncol = 2,
  nrow = 2)
