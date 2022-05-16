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
  geom_lineribbon(color = brewer_line_color, step = "hv", key_glyph = "rect") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = my_labeller["R₀_t"], limits = c(0.25, 3)) +
  ggtitle("Posterior Basic Reproduction Number") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  my_theme +
  theme(legend.position = c(0.1, 0.8), legend.direction="horizontal")

Rt_plot <- 
  vector_gq %>% 
  filter(name == "Rₜ_t") %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = brewer_line_color, step = "hv", key_glyph = "rect") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = my_labeller["Rₜ_t"], limits = c(0.25, 3)) +
  ggtitle("Posterior Effective Reproduction Number") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  my_theme +
  guides(fill = "none")


Rt_R0_plot <- plot_grid(R0_plot, Rt_plot, ncol = 1, align = "hv")
save_plot_target_asp(filename = "figures/advancement_slides/Rt_R0_plot.pdf", plot = Rt_R0_plot, base_asp = 16/9, base_height = 6)
