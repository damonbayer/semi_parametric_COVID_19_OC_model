# Shape of the pandemic
library(tidyverse)
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


vector_gq_path <- 
  tibble(full_path = dir_ls("results/tidy_vector_generated_quantities")) %>% 
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

vector_gq <- 
  read_csv(vector_gq_path) %>% 
  filter(date <= max_date)

deaths_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = vector_gq %>% 
                    filter(name == "D"),
                  mapping = aes(ymin = .lower, ymax = .upper),
                  color = brewer_line_color, key_glyph = "rect") +
  geom_line(data = dat_tidy %>% 
               filter(name == "cumulative_deaths")) +
  geom_point(data = dat_tidy %>% 
               filter(name == "cumulative_deaths")) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  ggtitle("Posterior Latent and Observed Deaths") +
  my_theme +
  theme(legend.position = c(0.1, 0.8))

cases_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = vector_gq %>% 
                    filter(name == "C"),
                  mapping = aes(ymin = .lower, ymax = .upper),
                  color = brewer_line_color, key_glyph = "rect") +
  geom_line(data = dat_tidy %>% 
              filter(name == "cumulative_cases")) +
  geom_point(data = dat_tidy %>% 
               filter(name == "cumulative_cases")) +
  scale_y_continuous(name = "Cases", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  ggtitle("Posterior Latent and Observed Cases") +
  my_theme +
  guides(fill = "none")

prevalence_plot <- 
  ggplot(mapping = aes(date, value)) +
  geom_lineribbon(data = vector_gq %>% 
                    filter(name == "prevalence"),
                  mapping = aes(ymin = .lower, ymax = .upper),
                  color = brewer_line_color, key_glyph = "rect") +
  scale_y_continuous(name = "Prevalence", labels = comma) +
  scale_x_date(name = "Date", date_breaks = date_breaks, date_labels = date_labels) +
  ggtitle("Posterior Latent Prevalence") +
  my_theme +
  guides(fill = "none")


shape_of_the_pandemic_plot <- plot_grid(deaths_plot, cases_plot, prevalence_plot, nrow = 1, align = "hv")

save_plot_target_asp(filename = "figures/advancement_slides/shape_of_the_pandemic_plot.pdf",
                     plot = shape_of_the_pandemic_plot,
                     base_asp = 1730/650,
                     base_height = 5)
