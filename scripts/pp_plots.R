library(gridExtra)
library(tidyverse)
library(tidybayes)
library(fs)
library(scales)
library(latex2exp)
library(cowplot)
source("src/plot_functions.R")

max_t <- 42

dat <- 
  read_csv("data/oc_data.csv") %>% 
  mutate(index = 1:n()) %>% 
  rename(date = end_date) %>% 
  select(-start_date) %>% 
  filter(time <= max_t)

max_date <- dat %>% filter(time == max_t) %>% pull(date)


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
  filter(max_t == 42) %>% 
  mutate(pp_data = map(full_path, read_csv)) %>% 
  select(-full_path) %>% 
  arrange(model_design) %>% 
  unnest(pp_data) %>% 
  filter(date <= max_date)

pp_test_pos_sensitivity_plot <- 
  all_pp %>% 
  filter(name == "test_positivity") %>% 
  filter(.width == 0.8) %>% 
  ggplot(aes(date, value)) +
  geom_lineribbon(mapping = aes(ymin = .lower,
                                ymax = .upper,
                                fill = as_factor(model_design)),
                  alpha = 0.1, show.legend = F) +
  geom_point(data = dat_tidy %>% filter(name == "test_positivity")) +
  scale_y_continuous(name = "Test Positivity", labels = percent) +
  scale_x_date(name = "Date") +
  ggtitle("Sensitivity Analysis - Posterior Predictive Test Positivity")

pp_death_sensitivity_plot <- 
  all_pp %>% 
  filter(name == "deaths") %>% 
  filter(.width == 0.8) %>% 
  ggplot(aes(date, value)) +
  geom_lineribbon(mapping = aes(ymin = .lower,
                                ymax = .upper,
                                fill = as_factor(model_design)),
                  alpha = 0.1, show.legend = F) +
  geom_point(data = dat_tidy %>% filter(name == "deaths")) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date") +
  ggtitle("Sensitivity Analysis - Posterior Predictive Deaths")



pp_test_pos_plot <-
  all_pp %>% 
    filter(constant_alpha == F,
           constant_IFR == F,
           constant_R0 == F,
           double_IFR_0 == F,
           half_alpha_0 == F,
           half_R0_0 == F,
           half_S_0 == F,
           use_seroprev == T,
           use_tests == T) %>% 
  filter(name == "test_positivity") %>% 
  ggplot(aes(date, value)) +
  geom_lineribbon(mapping = aes(ymin = .lower,
                                ymax = .upper)) +
  geom_point(data = dat_tidy %>% filter(name == "test_positivity")) +
  scale_y_continuous(name = "Test Positivity", labels = percent) +
  scale_x_date(name = "Date") +
  ggtitle("Posterior Predictive Test Positivity") +
  my_theme + theme(legend.position = "right")

pp_deaths_plot <-
  all_pp %>% 
  filter(constant_alpha == F,
         constant_IFR == F,
         constant_R0 == F,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F,
         use_seroprev == T,
         use_tests == T) %>% 
  filter(name == "deaths") %>% 
  ggplot(aes(date, value)) +
  geom_lineribbon(mapping = aes(ymin = .lower,
                                ymax = .upper)) +
  geom_point(data = dat_tidy %>% filter(name == "deaths")) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date") +
  ggtitle("Posterior Predictive Deaths") +
  my_theme + theme(legend.position = "right")


c("pp_death_sensitivity_plot", "pp_deaths_plot", "pp_test_pos_plot", 
  "pp_test_pos_sensitivity_plot") %>% 
  walk(~save_plot_target_asp(filename = path("figures/advancement_slides", ., ext = "pdf"),
                             plot = get(.),
                             base_asp = 32/9,
                             base_height = 4))
