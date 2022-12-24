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
  pivot_wider(names_from = name, values_from = value)

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
  filter(date <= max_date) %>% 
  filter(name %in% c("C", "I", "IFR_t", "R₀_t", "Rₜ_t", "α_t"))

R_y_lims <- c(0.25, 3.01)

R0_plot <- 
  all_posterior_generated_quantities %>% 
  filter(name == "R₀_t") %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = brewer_line_color, step = "hv", key_glyph = "rect") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = my_labeller["R₀_t"], limits = R_y_lims) +
  ggtitle("Posterior Basic Reproduction Number") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  my_theme +
  theme(legend.position = c(0.05, 0.8), legend.direction="vertical")

Rt_plot <- 
  all_posterior_generated_quantities %>% 
  filter(name == "Rₜ_t") %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = brewer_line_color, step = "hv", key_glyph = "rect") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = my_labeller["Rₜ_t"], limits = R_y_lims) +
  ggtitle("Posterior Effective Reproduction Number") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  my_theme +
  theme(legend.position = "none")

alpha_plot <- 
  all_posterior_generated_quantities %>% 
  filter(name == "α_t") %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = brewer_line_color, step = "hv", key_glyph = "rect") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = my_labeller["α_t"], limits = c(0, NA)) +
  ggtitle(TeX("Posterior $\\alpha$", bold = T)) +
  my_theme +
  theme(legend.position = "none")

ifr_t_plot <- 
  all_posterior_generated_quantities %>%
  filter(name == "IFR_t") %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = brewer_line_color, step = "hv", key_glyph = "rect") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = my_labeller["IFR_t"], limits = c(0, NA)) +
  ggtitle(TeX("Posterior $\\eta$", bold = T)) +
  my_theme +
  theme(legend.position = "none")


cumulative_latent_case_ratio_plot <- 
  all_posterior_generated_quantities %>% 
  filter(name == "C") %>% 
  left_join(dat_tidy) %>% 
  drop_na() %>% 
  mutate(across(c(value, .lower, .upper), ~`/`(., cumulative_cases))) %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = brewer_line_color, step = "hv", key_glyph = "rect") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = "Cumulative Latent:Case Ratio", limits = c(1, NA), breaks = seq(5, 30, 5)) +
  ggtitle("Posterior Cumulative Latent:Case Ratio") +
  my_theme +
  theme(legend.position = "none")

weekly_latent_case_ratio_plot <- 
  all_posterior_generated_quantities %>% 
  filter(name == "I") %>% 
  left_join(dat_tidy) %>% 
  drop_na() %>% 
  mutate(across(c(value, .lower, .upper), ~`/`(., weekly_cases))) %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(color = brewer_line_color, step = "hv", key_glyph = "rect") +
  scale_x_date(name = "Date", date_breaks = "3 months", date_labels = "%b %y") +
  scale_y_continuous(name = "Weekly Latent:Case Ratio", limits = c(1, NA), breaks = seq(5, 30, 5)) +
  ggtitle("Posterior Weekly Latent:Case Ratio") +
  my_theme +
  theme(legend.position = "none")


main_posterior_results_plot <- 
  plot_grid(R0_plot,
            Rt_plot,
            ifr_t_plot,
            alpha_plot,
            weekly_latent_case_ratio_plot,
            cumulative_latent_case_ratio_plot,
            ncol = 2, align = "hv")

save_plot(filename = path(figures_dir, "main_posterior_results_plot", ext = "pdf"),
          plot = main_posterior_results_plot,
          ncol = 2,
          nrow = 3,
          device = cairo_pdf
)
