library(tidyverse)
source("src/plot_functions.R")
library(latex2exp)
library(fs)

coverage_results <- read_csv("results/simulation/oc_like/coverage_results.csv")
shrinkage_results <- read_csv("results/simulation/oc_like/shrinkage_results.csv")

# Time-Varying Shrinkage --------------------------------------------------
generated_quantities_simulation_time_varying_shrinkage_plot <-
  shrinkage_results %>%
  filter(!is.na(time)) %>%
  filter(str_ends(name, "_t")) %>%
  ggplot(aes(time, shrinkage, color = name, ymin = .lower, ymax = .upper)) +
  geom_line() +
  geom_point() +
  geom_errorbar() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_discrete(name = "Parameter", labels = my_labeller) +
  scale_y_continuous(name = "Contraction", label = percent) +
  scale_x_continuous("Time") +
  ggtitle("Posterior Contraction Properties of Time-Varying Parameters",
    subtitle = "Median and 95% Interval from 200 simulations"
  ) +
  my_theme

# Compartment Shrinkage ---------------------------------------------------
generated_quantities_simulation_compartment_shrinkage_plot <-
  shrinkage_results %>%
  filter(!is.na(time)) %>%
  filter(name %in% c("S", "E", "I", "R", "D")) %>%
  mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))) %>%
  drop_na() %>%
  ggplot(aes(time, shrinkage, color = name, ymin = .lower, ymax = .upper)) +
  geom_line() +
  geom_point() +
  geom_errorbar() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_discrete(name = "Compartment", labels = my_labeller) +
  scale_y_continuous(name = "Contraction", label = percent) +
  scale_x_continuous("Time") +
  ggtitle("Posterior Contraction Properties of Compatment Sizes",
    subtitle = "Median and 95% Interval from 200 simulations"
  ) +
  my_theme


# Scalar Shrinkage --------------------------------------------------------
generated_quantities_simulation_scalar_shrinkage_plot <-
  shrinkage_results %>%
  filter(is.na(time)) %>%
  filter(!(name %in% c("S_SEI", "I_EI"))) %>%
  select(-time) %>%
  mutate(name = fct_reorder(name, shrinkage)) %>%
  ggplot(aes(shrinkage, name, xmin = .lower, xmax = .upper)) +
  geom_point() +
  geom_errorbar() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(name = "Parameter", labels = my_labeller) +
  scale_x_continuous(name = "Contraction", labels = percent) +
  ggtitle("Posterior Contraction Properties of Scalar Parameters",
          subtitle = "Median and 95% Interval from 200 simulations"
  ) +
  my_theme

# Time-Varying Coverage ---------------------------------------------------
generated_quantities_simulation_time_varying_coverage_plot <-
  coverage_results %>%
  filter(!is.na(time)) %>%
  filter(.width == 0.8) %>%
  filter(str_ends(name, "_t")) %>%
  ggplot(aes(time, mean, color = name, ymin = lower, ymax = upper)) +
  geom_line() +
  geom_point() +
  geom_errorbar() +
  scale_color_discrete(name = "Parameter", labels = my_labeller) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_y_continuous(name = "Coverage", label = percent) +
  scale_x_continuous("Time") +
  ggtitle("Posterior Coverage Properties of Time-Varying Parameters",
    subtitle = "Mean and 95% Confidence Interval from 200 Simulations"
  ) +
  my_theme

# Compartment Coverage ----------------------------------------------------
generated_quantities_simulation_compartment_coverage_plot <-
  coverage_results %>%
  filter(.width == 0.8) %>%
  filter(!is.na(time)) %>%
  filter(name %in% c("S", "E", "I", "R", "D")) %>%
  mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))) %>%
  ggplot(aes(time, mean, color = name, ymin = lower, ymax = upper)) +
  geom_line() +
  geom_point() +
  geom_errorbar() +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  scale_color_discrete(name = "Compartment", labels = my_labeller) +
  scale_y_continuous(name = "Coverage", label = percent) +
  scale_x_continuous("Time") +
  ggtitle("Posterior Coverage Properties of Compartment Sizes",
    subtitle = "Mean and 95% Confidence Interval from 200 Simulations"
  ) +
  my_theme

# Scalar Coverage ---------------------------------------------------------
generated_quantities_simulation_scalar_coverage_plot <-
  coverage_results %>%
  filter(.width == 0.8) %>%
  filter(is.na(time)) %>%
  filter(!(name %in% c("S_SEI", "I_EI"))) %>%
  select(-time) %>%
  mutate(name = fct_reorder(name, mean)) %>%
  ggplot(aes(mean, name, xmin = lower, xmax = upper)) +
  geom_point() +
  geom_errorbar() +
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  scale_y_discrete(name = "Parameter", labels = my_labeller) +
  scale_x_continuous(name = "Coverage", labels = percent) +
  ggtitle("Posterior Coverage Properties of Scalar Parameters",
          subtitle = "Mean and 95% Confidence Interval from 200 Simulations"
  ) +
  my_theme

# Save Plots --------------------------------------------------------------
c(
  "generated_quantities_simulation_compartment_coverage_plot",
  "generated_quantities_simulation_compartment_shrinkage_plot",
  "generated_quantities_simulation_time_varying_coverage_plot",
  "generated_quantities_simulation_time_varying_shrinkage_plot"
) |>
  walk(~ save_plot(
    filename = path(figures_dir, ., ext = "pdf"),
    plot = get(.),
    ncol = 1,
    nrow = 1, base_asp = 2
  ))

c(
  "generated_quantities_simulation_scalar_coverage_plot",
  "generated_quantities_simulation_scalar_shrinkage_plot"
) |>
  walk(~ save_plot(
    filename = path(figures_dir, ., ext = "pdf"),
    plot = get(.),
    base_asp = 1.75
  ))
