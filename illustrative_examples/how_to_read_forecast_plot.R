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

main_pp <- 
  all_pp %>% 
  filter(max_t == 42,
         constant_alpha == F,
         constant_IFR == F,
         constant_R0 == F,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F,
         use_seroprev == T,
         use_tests == T) %>% 
  unnest(pp_data) %>% 
  select(date, name, value, starts_with(".")) %>%
  mutate(forecast = date >= date %>% unique() %>% sort() %>% tail(4) %>% head(1))

first_forecast_date <- 
  main_pp %>% 
  filter(forecast) %>% 
  pull(date) %>% 
  min()

last_forecast_date <- 
  main_pp %>% 
  filter(forecast) %>% 
  pull(date) %>% 
  max()

last_modeling_date <- 
  main_pp %>% 
  filter(!forecast) %>% 
  pull(date) %>% 
  max()


first_modeling_date <- 
  main_pp %>% 
  filter(!forecast) %>% 
  pull(date) %>% 
  min()

how_to_read_forecast_forecast_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = main_pp %>% 
                    filter(name == "deaths"),
                  mapping = aes(ymin = .lower, ymax = .upper),
                  step = "mid",
                  color = brewer_line_color,
                  key_glyph = "rect") +
  geom_point(data = dat_tidy %>% 
               filter(name == "deaths") %>% 
               filter(date <= last_forecast_date)) +
  scale_x_date(name = "Date", limits = c(first_modeling_date, last_forecast_date)) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  my_theme +
  theme(legend.position = "none") +
  ggtitle("Posterior Deaths Forecast")

how_to_read_forecast_fit_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = main_pp %>% 
                    filter(name == "deaths") %>% 
                    filter(date <= last_modeling_date),
                  mapping = aes(ymin = .lower, ymax = .upper),
                  step = "mid",
                  key_glyph = "rect") +
  geom_point(data = dat_tidy %>% 
               filter(name == "deaths") %>% 
               filter(date <= last_forecast_date)) +
  scale_x_date(name = "Date", limits = c(first_modeling_date, last_forecast_date)) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  my_theme +
  scale_fill_brewer(palette = "Greens") +
  theme(legend.position = "none") +
  ggtitle("Posterior Deaths Forecast")


save_plot_target_asp(filename = "figures/advancement_slides/how_to_read_forecast_forecast_plot.pdf",
                     plot = how_to_read_forecast_forecast_plot,
                     base_asp = 16/9,
                     base_height = 5)

save_plot_target_asp(filename = "figures/advancement_slides/how_to_read_forecast_fit_plot.pdf",
                     plot = how_to_read_forecast_fit_plot,
                     base_asp = 16/9,
                     base_height = 5)
