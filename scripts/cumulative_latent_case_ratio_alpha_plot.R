# Case Ratio
# Reproductive Number
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
  unite(col = name, sum_type, name) %>% 
  filter(name == "cumulative_cases") %>% 
  select(date, cumulative_cases = value)

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
  filter(date <= max_date) %>% 
  filter(name %in% c("α_t", "C"))

cumulative_latent_case_ratio_plot <- 
  vector_gq %>% 
  filter(name == "C") %>% 
  left_join(dat_tidy) %>% 
  drop_na() %>% 
  mutate(across(c(value, .lower, .upper), ~`/`(., cumulative_cases))) %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = brewer_line_color, step = "hv") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = "Cumulative Latent:Case Ratio", limits = c(1, NA), breaks = seq(5, 30, 5)) +
  ggtitle("Posterior Cumulative Latent:Case Ratio") +
  my_theme +
  theme(legend.position = "right")

alpha_plot <- 
  vector_gq %>% 
  filter(name == "α_t") %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = brewer_line_color, step = "hv") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = my_labeller["α_t"], limits = c(0, NA)) +
  ggtitle(TeX("Posterior $\\alpha$", bold = T)) +
  my_theme +
  theme(legend.position = "right")


save_plot_target_asp(filename = "figures/advancement_slides/cumulative_latent_case_ratio_plot.pdf", plot = cumulative_latent_case_ratio_plot, base_asp = 32/9, base_height = 4)
save_plot_target_asp(filename = "figures/advancement_slides/alpha_plot.pdf", plot = alpha_plot, base_asp = 32/9, base_height = 4)
