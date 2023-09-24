library(tidyverse)
library(tidybayes)
library(scales)
source("src/plot_functions.R")

dat <-
  read_csv("illustrative_examples/age_structure/data/data.csv") %>%
  pivot_longer(-time)

true_generated_quantities <-
  read_csv("illustrative_examples/age_structure/data/true_generated_quantities.csv") %>%
  pivot_longer(-time, values_drop_na = TRUE) %>%
  add_count(name) %>%
  mutate(time = if_else(n == 1, NA_real_, time)) %>%
  select(-n) %>%
  mutate(
    population = if_else(
      str_ends(name, "_(o|y|total)"),
      str_extract(name, "[^_]+$"),
      NA_character_),
    name = if_else(
      str_ends(name, "_(o|y|total)"),
      str_extract(name, ".*(?=_)"),
      name
    )) %>%
  mutate(population = population %>%
           fct_relevel("y", "o", "total") %>%
           fct_recode(
             general = "y",
             vulnerable = "o",
             combined = "total"))

read_predictive <- function(predictive_file_path) {
  read_csv(predictive_file_path) %>%
    pivot_longer(-c(iteration, chain), names_to = "name_raw") %>%
    mutate(name =
             if_else(
               str_detect(name_raw, "\\[\\d+\\]"),
               str_extract(name_raw, ".+(?=\\[)"),
               name_raw),
           time = name_raw %>%
             str_extract("(?<=\\[)\\d+(?=\\])") %>%
             as.integer()) %>%
    mutate(name = str_remove(name, "new_")) %>%
    select(name, time, value) %>%
    group_by(name, time) %>%
    median_qi(.width = c(0.5, 0.8, 0.95)) %>%
    mutate(name = str_remove(name, "data_"))
}

read_generated_quantities <- function(generated_quantities_file_path) {
  read_csv(generated_quantities_file_path) %>%
    pivot_longer(-c(iteration, chain), names_to = "name_raw") %>%
    mutate(name =
             if_else(
               str_detect(name_raw, "\\[\\d+\\]"),
               str_extract(name_raw, ".+(?=\\[)"),
               name_raw),
           time = name_raw %>%
             str_extract("(?<=\\[)\\d+(?=\\])") %>%
             as.integer()) %>%
    mutate(time = if_else(name %in% c("cases_mean", "deaths_mean"), time, time - 1L)) %>%
    select(name, time, value) %>%
    group_by(name, time) %>%
    median_qi(.width = c(0.5, 0.8, 0.95))
}

prior_predictive <- read_predictive("illustrative_examples/age_structure/data/prior_predictive.csv")
posterior_predictive <- read_predictive("illustrative_examples/age_structure/data/posterior_predictive.csv")

prior_generated_quantities <- read_generated_quantities("illustrative_examples/age_structure/data/prior_generated_quantities.csv")
posterior_generated_quantities <- read_generated_quantities("illustrative_examples/age_structure/data/posterior_generated_quantities.csv")

all_generated_quantities <-
  bind_rows(mutate(prior_generated_quantities, distribution = "prior"),
            mutate(posterior_generated_quantities, distribution = "posterior")) %>%
  arrange(name, time, distribution)

true_generated_quantities %>%
  filter(!is.na(population)) %>%
  ggplot(aes(time, value, color = population)) +
  facet_wrap(~name, scales = "free_y") +
  geom_line() +
  scale_color_discrete(name = "Population", labels = str_to_title)

data_ifr_age_structure_plot <-
  ggplot(mapping = aes(time, value)) +
  facet_wrap(. ~ name, scale = "free_y", labeller = as_labeller(. %>% str_replace("_", " ") %>% str_to_title()), ncol = 1) +
  geom_line(mapping = aes(color = population),
            data = true_generated_quantities %>%
              filter(str_starts(name, "latent_")) %>%
              mutate(name = str_remove(name, "latent_"))) +
  geom_point(data = dat) +
  scale_y_continuous(name = "Count", labels = comma) +
  scale_x_continuous(name = "Time") +
  scale_color_discrete(name = "Population Group", label = str_to_title) +
  ggtitle("Latent and Observed Cases and Deaths") +
  theme(legend.position = "bottom")

posterior_predictive_ifr_age_structure_plot <-
  ggplot(mapping = aes(time, value)) +
  facet_wrap(. ~ name, scale = "free_y", labeller = as_labeller(. %>% str_replace("_", " ") %>% str_to_title()),
             nrow = 2) +
  geom_lineribbon(data = posterior_predictive,
                  mapping = aes(ymin = .lower, ymax = .upper),
                  step = "mid",
                  color = brewer_line_color) +
  geom_point(data = dat) +
  scale_y_continuous(name = "Count", labels = comma) +
  scale_x_continuous(name = "Time") +
  ggtitle("Posterior Predictive for Structured Population",
          "Modeled as Unstructured Population") +
  my_theme

ifr_age_structure_generated_quantities_simulation_compartment_plot <-
  ggplot(mapping = aes(time, value)) +
  facet_wrap(~name, scales = "free_y") +
  geom_lineribbon(mapping = aes(ymin = .lower, ymax =.upper, color = distribution, fill = distribution),
                  data = all_generated_quantities %>%
                    filter(.width == 0.8) %>%
                    filter(name %in% c("S", "E", "I", "R", "D")) %>%
                    mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))),
                  alpha = 0.5) +
  geom_line(data = true_generated_quantities %>%
              filter(name %in% c("S", "E", "I", "R", "D"),
                     population == "combined") %>%
              mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))),
            linetype = "dashed") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Count", labels = comma) +
  scale_color_discrete(name = "Distribution", label = str_to_title) +
  scale_fill_discrete(name = "Distribution", label = str_to_title) +
  ggtitle("Prior and Posterior Credible Intervals for Compartments",
          "Structured Population Modeled as Unstructured Population,\n80% credible intervals, true values in black") +
  theme(legend.position = c(5/6, 1/4))


ifr_age_structure_generated_quantities_simulation_time_varying_plot <-
  ggplot(mapping = aes(time, value)) +
  facet_wrap(~name, scales = "free_y", labeller = my_labeller_fn, ncol = 1) +
  geom_lineribbon(mapping = aes(ymin = .lower, ymax =.upper, color = distribution, fill = distribution),
                  data = all_generated_quantities %>%
                    filter(.width == 0.8) %>%
                    filter(name %in% c("IFR_t", "Rₜ_t")),
                  alpha = 0.5,
                  step = "hv") +
  geom_step(data = true_generated_quantities %>%
              filter(name %in% c("IFR_t", "Rₜ_t"),
                     time <= all_generated_quantities %>%
                       filter(.width == 0.8) %>%
                       filter(name %in% c("IFR_t", "Rₜ_t")) %>%
                       pull(time) %>%
                       max()),
            linetype = "dashed") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Value") +
  scale_color_discrete(name = "Distribution", label = str_to_title) +
  scale_fill_discrete(name = "Distribution", label = str_to_title) +
  ggtitle("Prior and Posterior Credible Intervals for Time-Varying Parameters",
          "Structured Population Modeled as Unstructured Population,\n80% credible intervals, true values in black") +
  theme(legend.position = "bottom")

ifr_age_structure_generated_quantities_simulation_scalar_plot <-
  ggplot() +
  facet_wrap(~name, scales = "free_x",
             labeller = my_labeller_fn) +
  geom_pointinterval(data = all_generated_quantities %>%
                       filter(name %in% c("dur_infectious_days", "dur_latent_days", "R₀", "β")),
                     mapping = aes(y = distribution, x = value, xmin = .lower, xmax = .upper, color = distribution)) +
  geom_vline(data = true_generated_quantities %>%
               filter(name %in% c("dur_infectious_days", "dur_latent_days", "R₀", "β")),
             mapping = aes(xintercept = value)) +
  scale_x_continuous(name = "Value") +
  scale_y_discrete(name = "Distribution", labels = str_to_title) +
  scale_fill_discrete(name = "Distribution") +
  scale_color_discrete(name = "Distribution") +
  ggtitle("Prior and Posterior Credible Intervals for Scalar Parameters",
          subtitle = "Structured Population Modeled as Unstructured Population,\n50%, 80%, 95% credible intervals, true values in black") +
  theme(legend.position = "none")


save_plot_target_asp(
  filename = path(defense_figures_dir, "data_ifr_age_structure_plot", ext = "pdf"),
  plot = data_ifr_age_structure_plot,
  ncol = 1,
  nrow = 2
)

save_plot_target_asp(
  filename = path(defense_figures_dir, "posterior_predictive_ifr_age_structure_plot", ext = "pdf"),
  plot = posterior_predictive_ifr_age_structure_plot,
  ncol = 1,
  nrow = 2
)

save_plot_target_asp(
  filename = path(defense_figures_dir, "ifr_age_structure_generated_quantities_simulation_scalar_plot", ext = "pdf"),
  plot = ifr_age_structure_generated_quantities_simulation_scalar_plot,
  ncol = 2,
  nrow = 2,
  base_asp = 2.2,
  base_height = 2
)

save_plot_target_asp(
  filename = path(defense_figures_dir, "ifr_age_structure_generated_quantities_simulation_time_varying_plot", ext = "pdf"),
  plot = ifr_age_structure_generated_quantities_simulation_time_varying_plot,
  ncol = 1,
  nrow = 2,
  base_asp = 2.2
)

save_plot_target_asp(
  filename = path(defense_figures_dir, "ifr_age_structure_generated_quantities_simulation_compartment_plot", ext = "pdf"),
  plot = ifr_age_structure_generated_quantities_simulation_compartment_plot,
  ncol = 3,
  nrow = 2,
  base_height = 2
)
