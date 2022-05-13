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
  filter(name %in% c("R₀_t", "Rₜ_t"))


R0_plot <- 
  vector_gq %>% 
  filter(name == "R₀_t") %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = brewer_line_color, step = "hv") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = my_labeller["R₀_t"], limits = c(0.25, 3)) +
  ggtitle("Posterior Basic Reproduction Number") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  my_theme +
  theme(legend.position = "right")

Rt_plot <- 
  vector_gq %>% 
  filter(name == "Rₜ_t") %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = brewer_line_color, step = "hv") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = my_labeller["Rₜ_tt"], limits = c(0.25, 3)) +
  ggtitle("Posterior Effective Reproduction Number") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  my_theme +
  theme(legend.position = "right")

save_plot_target_asp(filename = "figures/advancement_slides/R0_plot.pdf", plot = R0_plot, base_asp = 32/9, base_height = 4)
save_plot_target_asp(filename = "figures/advancement_slides/Rt_plot.pdf", plot = Rt_plot, base_asp = 32/9, base_height = 4)
