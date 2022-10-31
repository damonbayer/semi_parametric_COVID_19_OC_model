library(tidyverse)
source("src/plot_functions.R")

dat <- read_csv("data/oc_data.csv") %>% 
  filter(time <= 46)

date_breaks <- "3 months"
date_labels <- "%b %y"
line_size <- 1
point_size <- line_size + 1

binned_data_plot <-
  ggplot(dat, aes(end_date, tests)) +
  geom_line(size = line_size) +
  geom_point(size = point_size) +
  scale_y_continuous(name = "Tests", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  ggplot(dat, aes(end_date, cases)) +
  geom_line(size = line_size) +
  geom_point(size = point_size) +
  scale_y_continuous(name = "Cases", labels = comma) +
  scale_x_date(name = "Date", date_breaks, date_labels = date_labels) +
  ggplot(dat, aes(end_date, deaths)) +
  geom_line(size = line_size) +
  geom_point(size = point_size) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date", date_breaks, date_labels = date_labels) +
  ggplot(dat, aes(end_date, cases / tests)) +
  geom_line(size = line_size) +
  geom_point(size = point_size) +
  scale_y_continuous(name = "Testing Positivity",
                     labels = function(.) scales::percent(., accuracy = 1),
                     limits = c(0, NA)) +
  scale_x_date(name = "Date", date_breaks, date_labels = date_labels) +
  patchwork::plot_layout(ncol = 2, nrow = 2) +
  patchwork::plot_annotation(title = str_c("Orange County", ", CA Data"),
                             subtitle = "Counts binned into weekly periods")


save_plot(filename = path(figures_dir, "binned_data_plot", ext = "pdf"), plot = binned_data_plot, ncol = 2, nrow = 2)

simulated_dat <- read_csv("/Users/damon/Documents/semi_parametric_COVID_19_OC_model/data/simulated_data/simulated_data_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  slice(1) %>% 
  select(-c(iteration, chain, starts_with("data_seroprev_cases"))) %>% 
  pivot_longer(everything()) %>% 
  mutate(time = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  mutate(name = if_else(str_detect(name,"^.+\\["), str_extract(name, "^.+(?=\\[)"), name)) %>% 
  mutate(name = str_remove(name, "data_new_")) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  left_join(dat %>% select(-deaths, - cases))
