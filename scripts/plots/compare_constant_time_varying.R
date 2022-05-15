library(tidyverse)
library(tidybayes)
source("src/plot_functions.R")


library(gridExtra)
library(tidyverse)
library(tidybayes)
library(fs)
library(scales)
library(latex2exp)
library(cowplot)

max_t <- 42

dat <- 
  read_csv("data/oc_data.csv") %>% 
  mutate(index = 1:n()) %>% 
  rename(date = end_date) %>% 
  select(-start_date) %>% 
  filter(time <= max_t)

max_date <- dat %>% filter(time == max_t) %>% pull(date)

model_table <-
  read_csv("model_table.csv") %>% 
  select(-model_id, -seed) %>% 
  distinct()


all_vector_gq <- 
  tibble(full_path = dir_ls("results/tidy_vector_generated_quantities")) %>% 
  mutate(model_design = full_path %>% 
           str_extract("(?<=model_design=)\\d+") %>% 
           as.integer()) %>% 
  left_join(model_table) %>% 
  filter(max_t == 42) %>%
  filter(double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F,
         use_seroprev == T,
         use_tests == T) %>%
  select(full_path, model_design, starts_with("constant")) %>% 
  mutate("all_time_varying" = !constant_alpha & !constant_IFR & !constant_R0) %>% 
  pivot_longer(-c(full_path, model_design), names_to = "model_name") %>% 
  filter(value) %>% 
  select(-value) %>% 
  mutate(vector_gq_data = map(full_path, read_csv)) %>% 
  select(-full_path) %>% 
  arrange(model_design) %>% 
  unnest(vector_gq_data)

all_scalar_gq <- 
  tibble(full_path = dir_ls("results/tidy_scalar_generated_quantities")) %>% 
  mutate(model_design = full_path %>% 
           str_extract("(?<=model_design=)\\d+") %>% 
           as.integer()) %>% 
  left_join(model_table) %>% 
  filter(max_t == 42) %>%
  filter(double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F,
         use_seroprev == T,
         use_tests == T) %>%
  select(full_path, model_design, starts_with("constant")) %>% 
  mutate("all_time_varying" = !constant_alpha & !constant_IFR & !constant_R0) %>% 
  pivot_longer(-c(full_path, model_design), names_to = "model_name") %>% 
  filter(value) %>% 
  select(-value) %>% 
  mutate(scalar_gq_data = map(full_path, read_csv)) %>% 
  select(-full_path) %>% 
  arrange(model_design) %>% 
  unnest(scalar_gq_data)



all_gq <- 
  bind_rows(
  all_vector_gq %>% 
    filter(name %in% c("IFR_t", "R₀_t", "α_t")) %>% 
    mutate(name = str_remove(name, "_t")),
  all_scalar_gq %>% 
    filter(name %in% c("α", "R₀", "IFR")) %>% 
    right_join(distinct(all_vector_gq, model_design, date)) %>% 
    drop_na()
) %>% 
  filter(date <= max_date)

compare_constant_time_varying_plot <- 
  all_gq %>% 
  filter(.width == 0.8) %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper, fill = model_name)) +
  facet_wrap(. ~ name, scales = "free_y", labeller = my_labeller_fn) +
  geom_lineribbon(alpha = 0.5) +
  theme(legend.position = "bottom") +
  scale_y_continuous(name = "Value") +
  scale_x_date(name = "Date") +
  scale_fill_discrete(name = "Model",
                      labels = TeX(c(all_time_varying = "All Time Varying",
                                     constant_R0 = "Constant $R_0$",
                                     constant_alpha = "Constant $\\alpha",
                                     constant_IFR  = "Constant $\\eta$"))) +
  theme(legend.position = c(0.025, 3/4))

save_plot_target_asp(filename = path("figures/advancement_slides", "compare_constant_time_varying_plot", ext = "pdf"),
                           plot = compare_constant_time_varying_plot, ncol = 3, nrow = 1,
                           base_asp = 1730/650,
                     base_height = 5)
