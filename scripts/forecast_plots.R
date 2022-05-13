library(gridExtra)
library(tidyverse)
library(tidybayes)
library(fs)
library(scales)
library(latex2exp)
library(cowplot)
source("src/plot_functions.R")

dat <- 
  read_csv("data/oc_data.csv") %>% 
  mutate(index = 1:n()) %>% 
  rename(date = end_date) %>% 
  select(-start_date)

index_date_conversion <- 
  dat %>% 
  select(index, date)

dat_tidy <-
  dat %>% 
  select(date, cases, tests, deaths) %>% 
  mutate(test_positivity = cases / tests) %>% 
  select(-tests) %>% 
  pivot_longer(-date)

model_table <-
  read_csv("model_table.csv") %>% 
  select(-model_id, -seed) %>% 
  distinct()

all_pp <- 
  tibble(full_path = dir_ls("results/tidy_posterior_predictive")) %>% 
  mutate(model_design = full_path %>% 
           str_extract("(?<=model_design=)\\d+") %>% 
           as.integer()) %>% 
  left_join(model_table) %>% 
  mutate(pp_data = map(full_path, read_csv)) %>% 
  select(-full_path) %>% 
  arrange(model_design)

weeks_ahead_dat <- 
  all_pp %>% 
  left_join(index_date_conversion %>% 
              rename(max_t = index,
                     max_date = date)) %>% 
  unnest(pp_data) %>% 
  filter(date >= max_date) %>% 
  mutate(weeks_ahead = as.numeric((date - max_date) / 7))

weeks_ahead_dat_tmp <- 
  weeks_ahead_dat %>% 
  filter(constant_alpha == F,
         constant_IFR == F,
         constant_R0 == F,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F) %>% 
  filter(weeks_ahead %in% c(1, 4)) %>% 
  arrange(date) %>% 
  group_by(use_seroprev, use_tests) %>% 
  mutate(model_description = str_c("use_seroprev=", use_seroprev, "\n", "use_tests=", use_tests)) %>% 
  rename_with(.fn = ~str_replace_all(., "_", " ") %>% str_to_title(), .cols = c(use_seroprev, use_tests, weeks_ahead))

dat_tidy_tmp <- 
  crossing(dat_tidy,
         weeks_ahead_dat_tmp %>% 
           group_by(`Use Seroprev`,
                     `Use Tests`,
                    `Weeks Ahead`) %>% 
           summarize(min_date = min(date),
                     max_date = max(date))) %>% 
  filter(date >= min_date,
         date <= max_date) %>% 
  select(-ends_with("_date"))


test_pos_forecast_plot <- 
  weeks_ahead_dat_tmp %>% 
  filter(name == "test_positivity") %>% 
  ggplot(aes(date, value)) +
  facet_grid(`Weeks Ahead` ~ `Use Seroprev`+`Use Tests`, scale = "free", labeller = label_both) +
  geom_lineribbon(mapping = aes(ymin = .lower, ymax = .upper), step = "mid", color = brewer_line_color) +
  geom_point(data = dat_tidy_tmp %>% filter(name == "test_positivity")) +
  scale_y_continuous("Test Positivity", labels = percent) +
  scale_x_date("Date") +
  ggtitle("Test Positivity Forecasts") +
  my_theme

deaths_forecast_plot <- 
  weeks_ahead_dat_tmp %>% 
  filter(name == "deaths") %>% 
  ggplot(aes(date, value)) +
  facet_grid(`Weeks Ahead` ~ `Use Seroprev`+`Use Tests`, scale = "free", labeller = label_both) +
  geom_lineribbon(mapping = aes(ymin = .lower, ymax = .upper), step = "mid", color = brewer_line_color) +
  geom_point(data = dat_tidy_tmp %>% filter(name == "deaths")) +
  scale_y_continuous("Test Positivity", labels = comma) +
  scale_x_date("Date") +
  ggtitle("Deaths Forecasts") +
  my_theme

save_plot_target_asp(filename = "figures/advancement_slides/test_pos_forecast_plot.pdf", plot = test_pos_forecast_plot, ncol = 4, nrow = 2, base_asp = 16/9)
save_plot_target_asp(filename = "figures/advancement_slides/deaths_forecast_plot.pdf", plot = deaths_forecast_plot, ncol = 4, nrow = 2, base_asp = 16/9)
