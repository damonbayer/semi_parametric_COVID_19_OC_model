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
