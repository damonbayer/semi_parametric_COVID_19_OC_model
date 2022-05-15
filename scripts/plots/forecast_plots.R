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
  filter(weeks_ahead %in% 0:4) %>% 
  arrange(date) %>% 
  group_by(use_seroprev, use_tests) %>% 
  mutate(model_description = str_c("use_seroprev=", use_seroprev, "\n", "use_tests=", use_tests)) %>% 
  rename_with(.fn = ~str_replace_all(., "_", " ") %>% str_to_title(), .cols = c(use_seroprev, use_tests, weeks_ahead))

dat_tidy_tmp <- 
  crossing(dat_tidy,
         weeks_ahead_dat_tmp %>% 
           group_by(
             `Use Seroprev`,
             `Use Tests`,
             `Weeks Ahead`
             ) %>%
           summarize(min_date = min(date),
                     max_date = max(date))) %>% 
  filter(date >= min_date,
         date <= max_date) %>% 
  select(-ends_with("_date"))



make_weeks_ahead_forecast_comparison_plot <- function(weeks_ahead = 4) {
  cases_forecast_comparison_plot <- 
    weeks_ahead_dat_tmp %>% 
    filter(`Weeks Ahead` == weeks_ahead,
           `Use Seroprev` == T) %>% 
    filter(name == "cases") %>% 
    ggplot(aes(date, value)) +
    facet_wrap(. ~ `Use Tests`, labeller = label_both) +
    geom_lineribbon(mapping = aes(ymin = .lower, ymax = .upper), step = "mid", color = brewer_line_color, key_glyph = "rect") +
    geom_point(data = dat_tidy_tmp %>% filter(name == "cases", `Use Seroprev` == T, `Weeks Ahead` == weeks_ahead)) +
    scale_y_continuous("Cases", labels = comma) +
    scale_x_date("Date") +
    ggtitle(str_c(weeks_ahead, " Weeks Ahead Cases Forecasts")) +
    my_theme
  
  test_positivity_forecast_comparison_plot <- 
    weeks_ahead_dat_tmp %>% 
    filter(`Weeks Ahead` == weeks_ahead,
           `Use Seroprev` == T) %>% 
    filter(name == "test_positivity") %>% 
    ggplot(aes(date, value)) +
    facet_wrap(. ~ `Use Tests`, labeller = label_both) +
    geom_lineribbon(mapping = aes(ymin = .lower, ymax = .upper), step = "mid", color = brewer_line_color, key_glyph = "rect") +
    geom_point(data = dat_tidy_tmp %>% filter(name == "test_positivity", `Use Seroprev` == T, `Weeks Ahead` == weeks_ahead)) +
    scale_y_continuous("Test Positivity", labels = percent) +
    scale_x_date("Date") +
    ggtitle(str_c(weeks_ahead, " Weeks Ahead Test Positivity Forecasts")) +
    my_theme
  
  deaths_forecast_comparison_plot <- 
    weeks_ahead_dat_tmp %>% 
    filter(`Weeks Ahead` == weeks_ahead,
           `Use Seroprev` == T) %>% 
    filter(name == "deaths") %>% 
    ggplot(aes(date, value)) +
    facet_wrap(. ~ `Use Tests`, labeller = label_both) +
    geom_lineribbon(mapping = aes(ymin = .lower, ymax = .upper), step = "mid", color = brewer_line_color, key_glyph = "rect") +
    geom_point(data = dat_tidy_tmp %>% filter(name == "deaths", `Use Seroprev` == T, `Weeks Ahead` == weeks_ahead)) +
    scale_y_continuous("Deaths", labels = comma) +
    scale_x_date("Date") +
    ggtitle(str_c(weeks_ahead, " Weeks Ahead Deaths Forecasts")) +
    my_theme
  
  
  save_plot_target_asp(filename = str_c("figures/advancement_slides/", weeks_ahead, "_weeks_ahead_test_positivity_forecasts_comparison_plot.pdf"),
                       plot = test_positivity_forecast_comparison_plot,
                       base_asp = 16/9)
  save_plot_target_asp(filename = str_c("figures/advancement_slides/", weeks_ahead, "_weeks_ahead_cases_forecasts_comparison_plot.pdf"),
                       plot = cases_forecast_comparison_plot,
                       base_asp = 16/9)
  save_plot_target_asp(filename = str_c("figures/advancement_slides/", weeks_ahead, "_weeks_ahead_deaths_forecasts_comparison_plot.pdf"),
                       plot = deaths_forecast_comparison_plot,
                       base_asp = 16/9)

}

make_weeks_ahead_forecast_plot <- function(weeks_ahead = 4, use_tests = T) {
  dat_tidy_single_tmp <- 
    crossing(dat_tidy,
             weeks_ahead_dat_tmp %>% 
               filter(`Use Seroprev` == T,
                      `Use Tests` == T) %>% 
               group_by(`Weeks Ahead`) %>% 
               summarize(min_date = min(date),
                         max_date = max(date))) %>% 
    filter(date >= min_date,
           date <= max_date) %>% 
    select(-ends_with("_date"))
  
  deaths_forecast_single_plot <- 
    weeks_ahead_dat_tmp %>% 
    filter(name == "deaths",
           `Use Seroprev` == T,
           `Use Tests` == use_tests,
           `Weeks Ahead` == weeks_ahead) %>% 
    ggplot(aes(date, value)) +
    geom_lineribbon(mapping = aes(ymin = .lower, ymax = .upper), step = "mid", color = brewer_line_color, key_glyph = "rect") +
    geom_point(data = dat_tidy_single_tmp %>%
                 filter(name == "deaths", `Weeks Ahead` == weeks_ahead)) +
    scale_y_continuous("Deaths", labels = comma) +
    scale_x_date("Date") +
    ggtitle(str_c(weeks_ahead, " Weeks Ahead Deaths Forecasts"),
            subtitle = str_c(if_else(use_tests, "", "Not "), "using test counts")) +
    my_theme +
    theme(legend.position = c(0.05, 0.8))
  
  cases_forecast_single_plot <- 
    weeks_ahead_dat_tmp %>% 
    filter(name == "cases",
           `Use Seroprev` == T,
           `Use Tests` == use_tests,
           `Weeks Ahead` == weeks_ahead) %>% 
    ggplot(aes(date, value)) +
    geom_lineribbon(mapping = aes(ymin = .lower, ymax = .upper), step = "mid", color = brewer_line_color, key_glyph = "rect") +
    geom_point(data = dat_tidy_single_tmp %>%
                 filter(name == "cases", `Weeks Ahead` == weeks_ahead)) +
    scale_y_continuous("Cases", labels = comma) +
    scale_x_date("Date") +
    ggtitle(str_c(weeks_ahead, " Weeks Ahead Cases Forecasts"),
            subtitle = str_c(if_else(use_tests, "", "Not "), "using test counts")) +
    my_theme +
    theme(legend.position = c(0.05, 0.8))
  
  test_positivity_forecast_single_plot <- 
    weeks_ahead_dat_tmp %>% 
    filter(name == "test_positivity",
           `Use Seroprev` == T,
           `Use Tests` == use_tests,
           `Weeks Ahead` == weeks_ahead) %>% 
    ggplot(aes(date, value)) +
    geom_lineribbon(mapping = aes(ymin = .lower, ymax = .upper), step = "mid", color = brewer_line_color, key_glyph = "rect") +
    geom_point(data = dat_tidy_single_tmp %>%
                 filter(name == "test_positivity", `Weeks Ahead` == weeks_ahead)) +
    scale_y_continuous("Test Positivity", labels = percent) +
    scale_x_date("Date") +
    ggtitle(str_c(weeks_ahead, " Weeks Ahead Test Positivity Forecasts"),
            subtitle = str_c(if_else(use_tests, "", "Not "), "using test counts")) +
    my_theme +
    theme(legend.position = c(0.05, 0.8))
  
  save_plot_target_asp(filename = str_c("figures/advancement_slides/", weeks_ahead, "_weeks_ahead_cases_forecasts_", if_else(use_tests, "", "not_"), "using_test_counts_plot.pdf"),
                       plot = cases_forecast_single_plot,
                       base_asp = 16/9)
  
  save_plot_target_asp(filename = str_c("figures/advancement_slides/", weeks_ahead, "_weeks_ahead_deaths_forecasts_", if_else(use_tests, "", "not_"), "using_test_counts_plot.pdf"),
                       plot = deaths_forecast_single_plot,
                       base_asp = 16/9)
  
  save_plot_target_asp(filename = str_c("figures/advancement_slides/", weeks_ahead, "_weeks_ahead_test_positivity_forecasts_", if_else(use_tests, "", "not_"), "using_test_counts_plot.pdf"),
                       plot = test_positivity_forecast_single_plot,
                       base_asp = 16/9)
}
walk(4:0, make_weeks_ahead_forecast_plot)
walk(4:0, make_weeks_ahead_forecast_comparison_plot)
