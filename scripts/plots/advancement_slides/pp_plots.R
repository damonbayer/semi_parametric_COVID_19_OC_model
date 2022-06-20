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


main_pp <- 
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
  select(date, name, value, starts_with("."))


# Sensitivity Plots -------------------------------------------------------
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



c("pp_death_sensitivity_plot", "pp_test_pos_sensitivity_plot") %>% 
  walk(~save_plot_target_asp(filename = path("figures/advancement_slides", ., ext = "pdf"),
                             plot = get(.),
                             base_asp = 32/9,
                             base_height = 4))


# Abstract Plot -----------------------------------------------------------
tmp_xlim <- c(lubridate::ymd("2020-06-07"), lubridate::ymd("2020-08-23"))
tmp_ylim <- c(1402, 8287)
tmp_asp <- as.numeric((tmp_xlim[2] - tmp_xlim[1]) / (tmp_ylim[2] - tmp_ylim[1]))


pp_abstract_plot <- 
  all_pp %>% 
  filter(name == "cases") %>% 
  filter(model_design == 69) %>%
  ggplot(aes(date, value)) +
  geom_lineribbon(mapping = aes(ymin = .lower, ymax = .upper),
                  color = brewer_line_color) +
  geom_point(data = dat_tidy %>%
               filter(name == "cases")) +
  my_theme +
  theme_nothing() +
  coord_fixed(ratio = tmp_asp, xlim = tmp_xlim, ylim = tmp_ylim)

save_plot(filename = "~/Desktop/pp_abstract.pdf", plot = pp_abstract.pdf, base_asp = 1)


# Main Plot ---------------------------------------------------------------
library(patchwork)
library(tidyverse)
source("src/plot_functions.R")

date_breaks <- "3 months"
date_labels <- "%b %y"
line_size <- 1
point_size <- line_size 

pp_tests_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_line(data = dat_tidy %>% 
              filter(name == "tests"))+
  geom_point(data = dat_tidy %>% 
               filter(name == "tests")) +
  scale_y_continuous(name = "Tests", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  guides(fill = "none")


pp_cases_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = main_pp %>% 
                    filter(name == "cases"),
                  mapping = aes(ymin = .lower,
                                ymax = .upper),
                  color = brewer_line_color, key_glyph = "rect") +
  geom_point(data = dat_tidy %>% 
               filter(name == "cases")) +
  scale_y_continuous(name = "Cases", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  my_theme +
  theme(legend.position = c(0.25, 0.75))

pp_deaths_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = main_pp %>% 
                    filter(name == "deaths"),
                  mapping = aes(ymin = .lower,
                                ymax = .upper),
                  color = brewer_line_color) +
  geom_point(data = dat_tidy %>% 
               filter(name == "deaths")) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  scale_fill_brewer() +
  guides(fill = "none")

pp_test_pos_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = main_pp %>% 
                    filter(name == "test_positivity"),
                  mapping = aes(ymin = .lower,
                                ymax = .upper),
                  color = brewer_line_color,) +
  geom_point(data = dat_tidy %>% 
               filter(name == "test_positivity")) +
  scale_y_continuous(name = "Testing Positivity", labels = percent) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  scale_fill_brewer() +
  guides(fill = "none")

main_pp_plot <- 
  pp_tests_plot +
  pp_cases_plot +
  pp_deaths_plot +
  pp_test_pos_plot +
  patchwork::plot_layout(ncol = 2, nrow = 2) +
  patchwork::plot_annotation(title = "Posterior Predictive")

save_plot_target_asp(filename = "figures/advancement_slides/main_pp_plot.pdf",
                     plot = main_pp_plot,
                     ncol = 2,
                     nrow = 2,
                     base_asp = 16/9, base_height = 3)
