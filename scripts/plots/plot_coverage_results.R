library(tidyverse)
source("src/plot_functions.R")
library(latex2exp)
library(fs)

generated_quantities_coverage_summary <- read_csv("results/generated_quantities_coverage_summary.csv")


# Time-Varying ------------------------------------------------------------
gq_simulation_time_varying_shrinkage_plot <- 
  generated_quantities_coverage_summary %>% 
  select(-ends_with("_error")) %>% 
  filter(!is.na(time)) %>% 
  filter(str_ends(name, "_t")) %>% 
  ggplot(aes(time, shrinkage, color = name)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0,  linetype = "dashed") +
  scale_color_discrete(name = "Parameter", labels = my_labeller) +
  scale_y_continuous(name = "Shrinkage", label = percent) +
  scale_x_continuous("Time") +
  ggtitle("Posterior Shrinkage Properties of Time-Varying Parameters",
          subtitle = "from 200 simulations") +
  my_theme

gq_simulation_time_varying_coverage_plot <- 
  generated_quantities_coverage_summary %>% 
  select(-ends_with("_error")) %>% 
  filter(!is.na(time)) %>% 
  filter(str_ends(name, "_t")) %>% 
  ggplot(aes(time, covered, color = name)) +
  scale_color_discrete(name = "Parameter", labels = my_labeller) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_y_continuous(name = "Coverage", label = percent) +
  scale_x_continuous("Time") +
  ggtitle("Posterior Coverage Properties of Time-Varying Parameters",
          subtitle = "from 200 simulations") +
  my_theme

# Shrinkage ---------------------------------------------------------------
gq_simulation_compartment_shrinkage_plot <- 
  generated_quantities_coverage_summary %>% 
  select(-ends_with("_error")) %>% 
  filter(!is.na(time)) %>% 
  filter(name %in% c("S", "E", "I", "R", "D")) %>% 
  mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))) %>% 
  ggplot(aes(time, shrinkage, color = name)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0,  linetype = "dashed")
  scale_color_discrete(name = "Compartment", labels = my_labeller) +
  scale_y_continuous(name = "Shrinkage", label = percent) +
  scale_x_continuous("Time") +
  ggtitle("Posterior Shrinkage Properties of Compatment Sizes",
          subtitle = "from 200 simulations") +
  my_theme

gq_simulation_compartment_coverage_plot <- 
  generated_quantities_coverage_summary %>% 
  select(-ends_with("_error")) %>% 
  filter(!is.na(time)) %>% 
  filter(name %in% c("S", "E", "I", "R", "D")) %>% 
  mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))) %>% 
  ggplot(aes(time, covered, color = name)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_color_discrete(name = "Compartment", labels = my_labeller) +
  scale_y_continuous(name = "Coverage", label = percent) +
  scale_x_continuous("Time") +
  ggtitle("Posterior Coverage Properties of Compartment Sizes",
          subtitle = "from 200 simulations") +
  my_theme

# Scalar ------------------------------------------------------------------

gq_simulation_scalar_shrinkage_plot <- 
  generated_quantities_coverage_summary %>% 
  select(-ends_with("_error")) %>% 
  filter(is.na(time)) %>% 
  filter(!(name %in%  c("S_SEI", "I_EI"))) %>% 
  select(-time) %>% 
  mutate(name = fct_reorder(name, shrinkage)) %>% 
  ggplot(aes(shrinkage, name)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(name = "Parameter", labels = my_labeller) +
  scale_x_continuous(name = "Shrinkage", labels = percent)  +
  ggtitle("Posterior Shrinkage Properties",
          subtitle = "Time-stationary parameters, 200 simulations") +
  my_theme

gq_simulation_scalar_coverage_plot <- 
  generated_quantities_coverage_summary %>% 
  select(-ends_with("_error")) %>% 
  filter(is.na(time)) %>% 
  filter(!(name %in%  c("S_SEI", "I_EI"))) %>% 
  select(-time) %>% 
  mutate(name = fct_reorder(name, covered)) %>% 
  ggplot(aes(covered, name)) +
  geom_point() +
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  scale_y_discrete(name = "Parameter", labels = my_labeller) +
  scale_x_continuous(name = "Coverage", labels = percent) +
  ggtitle("Posterior Coverage Properties",
          subtitle = "Time-stationary parameters, 200 simulations") +
  my_theme


# Save Plots --------------------------------------------------------------
c("gq_simulation_compartment_coverage_plot", "gq_simulation_compartment_shrinkage_plot", 
  "gq_simulation_time_varying_coverage_plot", "gq_simulation_time_varying_shrinkage_plot"
) %>% 
  walk(~save_plot_target_asp(filename = path("figures/advancement_slides", ., ext = "pdf"),
                             plot = get(.),
                             base_height = 3,
                             base_asp = 24/9))


c("gq_simulation_scalar_coverage_plot", "gq_simulation_scalar_shrinkage_plot") %>% 
  walk(~save_plot_target_asp(filename = path("figures/advancement_slides", ., ext = "pdf"),
                             plot = get(.),
                             base_height = 5,
                             base_asp = 8/9))

