library(tidyverse)
source("src/plot_functions.R")
model_table <-
  read_csv("model_table.csv") %>% 
  select(-model_id, -seed) %>% 
  distinct()

all_vector_gq <- 
  bind_rows(
    tibble(full_path = dir_ls("results/tidy_vector_generated_quantities"),
           model_type = "posterior"),
    tibble(full_path = dir_ls("results/tidy_vector_prior_generated_quantities"),
           model_type = "prior")) %>% 
  mutate(model_design = full_path %>% 
           str_extract("(?<=model_design=)\\d+") %>% 
           as.integer()) %>% 
  left_join(model_table) %>% 
  filter(max_t == max(max_t),
         !constant_alpha,
         !constant_IFR,
         !constant_R0,
         use_tests,
         use_seroprev) %>%
  select(full_path, model_type, model_design, starts_with("half"), starts_with("double")) %>% 
  rowwise() %>% 
  mutate(original = !any(c(half_alpha_0, half_R0_0, half_S_0, double_IFR_0))) %>% 
  pivot_longer(c(half_alpha_0, half_R0_0, half_S_0, double_IFR_0, original)) %>% 
  group_by(full_path, model_type, model_design) %>% 
  summarize(model_name = name[value], .groups = "drop") %>% 
  mutate(scalar_gq_data = map(full_path, read_csv)) %>% 
  select(-full_path) %>% 
  arrange(model_type, model_design) %>% 
  unnest(scalar_gq_data)


all_scalar_gq <-
  bind_rows(
    tibble(full_path = dir_ls("results/tidy_scalar_generated_quantities"),
           model_type = "posterior"),
    tibble(full_path = dir_ls("results/tidy_scalar_prior_generated_quantities"),
           model_type = "prior")) %>% 
  mutate(model_design = full_path %>% 
           str_extract("(?<=model_design=)\\d+") %>% 
           as.integer()) %>% 
  left_join(model_table) %>% 
  filter(max_t == max(max_t),
         !constant_alpha,
         !constant_IFR,
         !constant_R0,
         use_tests,
         use_seroprev) %>%
  select(full_path, model_type, model_design, starts_with("half"), starts_with("double")) %>% 
  rowwise() %>% 
  mutate(original = !any(c(half_alpha_0, half_R0_0, half_S_0, double_IFR_0))) %>% 
  pivot_longer(c(half_alpha_0, half_R0_0, half_S_0, double_IFR_0, original)) %>% 
  group_by(full_path, model_type, model_design) %>% 
  summarize(model_name = name[value], .groups = "drop") %>% 
  mutate(scalar_gq_data = map(full_path, read_csv)) %>% 
  select(-full_path) %>% 
  arrange(model_type, model_design) %>% 
  unnest(scalar_gq_data)


model_name_converter <- 
  tibble(
    model_name = c("original", "half_alpha_0", "half_R0_0", "half_S_0", "double_IFR_0"),
    model_name_to_plot = c("Original", "Half $\\exp(\\tilde{\\alpha}_1)$", "Half $\\exp(\\tilde{R}_{0,1})$", "Half $S_0$", "Double $\\expit(\\tilde{\\eta}_1)$"))

# all_scalar_gq %>% 
#   left_join(model_name_converter) %>% 
#   mutate(name_for_plot = str_c(model_name_to_plot, " (", str_to_title(model_type), ")")) %>% 
#   ggplot(aes(value, name_for_plot, xmin = .lower, xmax = .upper, color = model_type)) +
#   facet_wrap(. ~ name, scales = "free_x", labeller = my_labeller_fn) +
#   geom_interval(alpha = 1/2) +
#   scale_y_discrete(labels = ~TeX(.)) +
#   scale_color_discrete(name = "Distribution", label = str_to_title) +
#   theme(legend.position = c(3/5, 1/4))

scalar_sensitivity_plot <- 
  all_scalar_gq %>% 
  left_join(model_name_converter) %>% 
  filter(.width == 0.8) %>%
  ggplot(aes(value, model_name_to_plot, xmin = .lower, xmax = .upper, color = model_type)) +
  facet_wrap(. ~ name, scales = "free_x", labeller = my_labeller_fn) +
  geom_interval() * 
  blend("lighten") +
  geom_interval() * 
  blend("multiply", alpha = 0.6) +
  scale_y_discrete("Model", labels = ~TeX(.)) +
  scale_color_discrete(name = "Distribution", label = str_to_title) +
  theme(legend.position = c(3/5, 1/4)) +
  scale_x_continuous("Value") +
  ggtitle("Sensitivity Analysis Results",
          subtitle = "80% Credible Intervals")

compartments_sensitivity_plot <- 
  all_vector_gq %>% 
  left_join(model_name_converter) %>% 
  filter(name %in% c("S", "E", "I", "R", "D")) %>% 
  filter(.width == 0.8) %>% 
  mutate(name = name %>% fct_relevel(c("S", "E", "I", "R", "D"))) %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper, fill = model_type)) +
  facet_grid(name ~ model_name_to_plot, scales = "free_y", labeller = as_labeller(function(string) TeX(string), label_parsed)) +
  geom_lineribbon() *
  blend("lighten") +
  geom_lineribbon() *
  blend("multiply", alpha = 0.6) +
  scale_fill_discrete(name = "Distribution", label = str_to_title) +
  scale_y_continuous(name = "Count", labels = comma) +
  scale_x_date("Date") +
  theme(legend.position = "bottom") +
  ggtitle("Sensitivity Analysis Results",
          subtitle = "80% Credible Intervals")
 
time_varying_sensitivity_plot <- 
  all_vector_gq %>% 
  left_join(model_name_converter) %>% 
  filter(name %in% c("IFR_t", "R₀_t", "Rₜ_t", "α_t", "β_t")) %>%
  filter(.width == 0.8) %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper, fill = model_type)) +
  facet_grid(name ~ model_name_to_plot, scales = "free_y", labeller = labeller(name = my_labeller_fn, model_name_to_plot = as_labeller(function(string) TeX(string), label_parsed))) +
  geom_lineribbon() *
  blend("lighten") +
  geom_lineribbon() *
  blend("multiply", alpha = 0.6) +
  scale_fill_discrete(name = "Distribution", label = str_to_title) +
  scale_y_continuous(name = "Value", labels = comma) +
  scale_x_date("Date") +
  theme(legend.position = "bottom") +
  ggtitle("Sensitivity Analysis Results",
        subtitle = "80% Credible Intervals")

save_plot(filename = path(figures_dir, "scalar_sensitivity_plot", ext = "pdf"),
          plot = scalar_sensitivity_plot,
          ncol = 4, nrow = 3, base_asp = 1,
          device = cairo_pdf)

save_plot(filename = path(figures_dir, "compartments_sensitivity_plot", ext = "pdf"),
          plot = compartments_sensitivity_plot,
          ncol = 5, nrow = 5,
          device = cairo_pdf)

save_plot(filename = path(figures_dir, "time_varying_sensitivity_plot", ext = "pdf"),
          plot = time_varying_sensitivity_plot,
          ncol = 5, nrow = 5,
          device = cairo_pdf)
