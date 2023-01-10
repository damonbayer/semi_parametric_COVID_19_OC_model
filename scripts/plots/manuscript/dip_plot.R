library(tidyverse)
library(lubridate)
library(glue)
source("src/plot_functions.R")
model_table <- read_csv("model_table.csv")

max_t <- 42

date_breaks = "3 months"
date_labels = "%b %y"

dat <- read_csv("data/oc_data.csv") %>% 
  filter(time <= max_t)

max_date <- dat %>% filter(time == max_t) %>% pull(end_date)

dat_tidy <- 
  dat %>% 
  select(date = end_date, cases, tests, deaths) %>% 
  pivot_longer(-date, values_to = "weekly") %>% 
  group_by(name) %>% 
  mutate(cumulative = cumsum(weekly)) %>% 
  pivot_longer(cols = c(weekly, cumulative), names_to = "sum_type") %>% 
  unite(col = name, sum_type, name)


posterior_generated_quantities_path <- 
  tibble(full_path = dir_ls("results/tidy_posterior_generated_quantities")) %>% 
  mutate(model_id = full_path %>% 
           str_extract("(?<=model_id=)\\d+") %>% 
           as.numeric()) %>% 
  left_join(model_table %>% distinct(model_design, .keep_all = T)) %>% 
  select(-model_id, -seed) %>% 
  filter(constant_alpha == F,
         constant_IFR == F,
         constant_R0 == F,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F,
         max_t == 42,
         use_seroprev == T,
         use_tests == T) %>% 
  pull(full_path)

all_posterior_generated_quantities <- 
  read_csv(posterior_generated_quantities_path) %>% 
  filter(date <= max_date)

deaths_annotation_data_x <- lubridate::ymd("2020-12-31")

deaths_annotation_data_x_end <- max(dat_tidy$date) - 7

deaths_annotation_data_y_end <- 
  dat_tidy %>%
  filter(date == deaths_annotation_data_x_end,
         name == "cumulative_deaths") %>%
  pull(value)

deaths_annotation_data_y <- 500

deaths_annotation_posterior_x <- lubridate::ymd("2020-06-01")

deaths_annotation_posterior_x_end <- max(all_posterior_generated_quantities$date) - 7

deaths_annotation_posterior_y_end <- 
  all_posterior_generated_quantities %>% 
  filter(name == "D",
         date == deaths_annotation_posterior_x_end,
         .width == 0.95) %>% 
  pull(.upper)

deaths_annotation_posterior_y <- 1000

deaths_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = all_posterior_generated_quantities %>% 
                    filter(name == "D"),
                  mapping = aes(ymin = .lower, ymax = .upper),
                  color = brewer_line_color, key_glyph = "rect") +
  geom_line(data = dat_tidy %>% 
              filter(name == "cumulative_deaths")) +
  geom_point(data = dat_tidy %>% 
               filter(name == "cumulative_deaths")) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  annotate(geom = "curve",
           x = deaths_annotation_data_x,
           xend = deaths_annotation_data_x_end + 3,
           y = deaths_annotation_data_y,
           yend = deaths_annotation_data_y_end * 0.98,
           curvature = 0.25,
           arrow = arrow(length = unit(4, "mm"))) +
  annotate(geom = "text",
           x = deaths_annotation_data_x,
           y = deaths_annotation_data_y,
           label = "Observed Data",
           hjust = "right",
           vjust = "center") +
  annotate(geom = "curve",
           x = deaths_annotation_posterior_x,
           xend = deaths_annotation_posterior_x_end,
           y = deaths_annotation_posterior_y,
           yend = deaths_annotation_posterior_y_end,
           curvature = -0.25,
           arrow = arrow(length = unit(4, "mm"))) +
  annotate(geom = "text",
           x = deaths_annotation_posterior_x,
           y = deaths_annotation_posterior_y,
           label = "Posterior",
           hjust = "center",
           vjust = "top") +
  ggtitle("Posterior Latent and Observed Deaths") +
  my_theme +
  theme(legend.position = c(0.1, 0.8))

cases_annotation_x <- lubridate::ymd("2020-07-01")

cases_annotation_x_end <- lubridate::ymd("2020-08-16")

cases_annotation_y <- 1250000

cases_annotation_y_end <- 
  all_posterior_generated_quantities %>% 
  filter(name == "C",
         date == cases_annotation_x_end,
         .width == 0.95) %>% 
  pull(.upper)

seroprev_interval <-
  all_posterior_generated_quantities %>% 
  filter(name == "C",
         date == cases_annotation_x_end,
         .width == 0.95) %>% 
  select(value, .lower, .upper) %>% 
  unlist() %>% 
  `/`(popsize) %>% 
  percent()

cases_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = all_posterior_generated_quantities %>% 
                    filter(name == "C"),
                  mapping = aes(ymin = .lower, ymax = .upper),
                  color = brewer_line_color, key_glyph = "rect") +
  geom_line(data = dat_tidy %>% 
              filter(name == "cumulative_cases")) +
  geom_point(data = dat_tidy %>% 
               filter(name == "cumulative_cases")) +
  annotate(geom = "curve",
           x = cases_annotation_x,
           xend = cases_annotation_x_end,
           y = cases_annotation_y,
           yend = cases_annotation_y_end,
           curvature = 0.25,
           arrow = arrow(length = unit(4, "mm"))) +
  annotate(geom = "text",
           x = cases_annotation_x,
           y = cases_annotation_y * 1.01,
           label = glue("Posterior\nSeroprevalence:\n{seroprev_interval[['value']]}\n({seroprev_interval[['.lower']]}, {seroprev_interval[['.upper']]})"),
           hjust = "center",
           vjust = "bottom") +
  scale_y_continuous(name = "Cases", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  ggtitle("Posterior Latent and Observed Cases") +
  my_theme +
  guides(fill = "none")

prevalence_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = all_posterior_generated_quantities %>% 
                    filter(name == "prevalence"),
                  mapping = aes(ymin = .lower, ymax = .upper),
                  color = brewer_line_color, key_glyph = "rect") +
  scale_y_continuous(name = "Prevalence", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  ggtitle("Posterior Latent Prevalence") +
  my_theme +
  guides(fill = "none")


dip_plot <- plot_grid(deaths_plot, cases_plot, prevalence_plot, nrow = 1, align = "hv")

save_plot(filename = path(figures_dir, "dip_plot", ext = "pdf"),
          plot = dip_plot,
          ncol = 3,
          nrow = 1,
          base_asp = 1.25, base_height = 3.5)
