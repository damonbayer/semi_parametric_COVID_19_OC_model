library(tidyverse)
library(tidybayes)
library(fs)
source("src/plot_functions.R")

target_sim_id <- 1

real_dat <- read_csv("data/oc_data.csv") %>% 
  filter(time <= 42)

simulated_dat <- 
  read_csv("data/simulation/oc_like/simulated_data.csv") %>% 
  slice(target_sim_id) %>% 
  select(-c(iteration, chain, starts_with("data_seroprev_cases"))) %>% 
  pivot_longer(everything()) %>% 
  mutate(time = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  mutate(name = if_else(str_detect(name,"^.+\\["), str_extract(name, "^.+(?=\\[)"), name)) %>% 
  mutate(name = str_remove(name, "data_new_")) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  left_join(real_dat %>% select(-deaths, - cases))

line_size <- 1
point_size <- line_size + 1

simulated_binned_data_plot <-
  ggplot(simulated_dat, aes(time, tests)) +
  geom_line(size = line_size) +
  geom_point(size = point_size) +
  scale_y_continuous(name = "Tests", labels = comma) +
  scale_x_continuous(name = "Time") +
  ggplot(simulated_dat, aes(time, cases)) +
  geom_line(size = line_size) +
  geom_point(size = point_size) +
  scale_y_continuous(name = "Cases", labels = comma) +
  scale_x_continuous(name = "Time") +
  ggplot(simulated_dat, aes(time, deaths)) +
  geom_line(size = line_size) +
  geom_point(size = point_size) +
  scale_y_continuous(name = "Deaths", labels = comma) +
  scale_x_continuous(name = "Time") +
  ggplot(simulated_dat, aes(time, cases / tests)) +
  geom_line(size = line_size) +
  geom_point(size = point_size) +
  scale_y_continuous(name = "Testing Positivity",
                     labels = function(.) scales::percent(., accuracy = 1),
                     limits = c(0, NA)) +
  scale_x_continuous(name = "Time") +
  patchwork::plot_layout(ncol = 2, nrow = 2) +
  patchwork::plot_annotation(title = str_c("Simulated data"),
                             subtitle = "Counts binned into weekly periods")

line_size <- 2

true_generated_quantities <-
  read_csv("data/simulation/oc_like/true_generated_quantities.csv") %>% 
  select(-iteration, -chain) %>% 
  pivot_longer(everything(), values_to = "true_value") %>% 
  mutate(time = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  mutate(name = if_else(str_detect(name,"^.+\\["), str_extract(name, "^.+(?=\\[)"), name)) %>% 
  mutate(time = time - 1)

tidy_prior_generated_quantities <- read_csv("results/simulation/oc_like/tidy_prior_generated_quantities/tidy_prior_generated_quantities.csv")

tidy_posterior_generated_quantities <- 
  tibble(file_path = dir_ls("results/simulation/oc_like/tidy_posterior_generated_quantities")) %>% 
  mutate(sim_id = file_path %>% str_extract("(?<=sim_id=)\\d+") %>% as.integer()) %>% 
  filter(sim_id == target_sim_id) %>% 
  pull(file_path) %>% 
  read_csv()

all_tidy_generated_quantities <- 
  bind_rows(
    mutate(tidy_prior_generated_quantities, distribution = "prior"),
    mutate(tidy_posterior_generated_quantities, distribution = "posterior"))


single_generated_quantities_simulation_scalar_plot <-
  ggplot() +
  geom_pointinterval(data = all_tidy_generated_quantities %>%
                       filter(is.na(time)) %>% 
                       filter(!(name %in%  c("S_SEI", "I_EI"))),
                     mapping = aes(y = distribution, x = value, xmin = .lower, xmax = .upper, color = distribution)) +
  geom_vline(data = true_generated_quantities %>%
               filter(is.na(time)) %>% 
               filter(!(name %in%  c("S_SEI", "I_EI"))),
             mapping = aes(xintercept = true_value),
             size = line_size) +
  facet_wrap(~name, scales = "free_x",
             labeller = my_labeller_fn) +
  scale_x_continuous(name = "Value") +
  scale_y_discrete(name = "Distribution", labels = str_to_title) +
  scale_fill_discrete(name = "Distribution") +
  scale_color_discrete(name = "Distribution") +
  ggtitle("Prior and Posterior Credible Intervals for Scalar Parameters",
          subtitle = "One simulated dataset, 50%, 80%, 95% credible intervals, true values in black") +
  theme(legend.position = "none")


single_generated_quantities_simulation_compartment_plot <-
  ggplot() +
  geom_lineribbon(data = all_tidy_generated_quantities %>%
                    filter(.width == 0.8) %>% 
                    filter(!is.na(time)) %>% 
                    filter(name %in% c("S", "E", "I", "R", "D")) %>% 
                    mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))),
                  mapping = aes(x= time, y = value, ymin = .lower, ymax = .upper, fill = distribution, color = distribution, group = .width),
                  alpha = 0.5,
                  key_glyph = "rect") +
  geom_point(data = true_generated_quantities %>%
               filter(!is.na(time)) %>%
               filter(name %in% c("S", "E", "I", "R", "D")) %>%
               mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))),
             mapping = aes(x = time, y = true_value)) +
  geom_line(data = true_generated_quantities %>%
              filter(!is.na(time)) %>%
              filter(name %in% c("S", "E", "I", "R", "D")) %>%
              mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))),
            mapping = aes(x = time, y = true_value)) +
  facet_wrap(. ~ name, scales = "free_y",
             labeller = my_labeller_fn) +
  scale_y_continuous(name = "Count", labels = comma) +
  scale_x_continuous(name = "Time") +
  scale_fill_discrete(name = "Distribution", labels = str_to_title) +
  scale_color_discrete(name = "Distribution", labels = str_to_title) +
  ggtitle("Prior and Posterior Credible Intervals for Compartments",
          subtitle = "One simulated dataset, 80% credible intervals, true values in black") +
  theme(legend.position = c(5/6, 1/4))

single_generated_quantities_simulation_time_varying_plot <-
  ggplot() +
  geom_lineribbon(data = all_tidy_generated_quantities %>%
                    filter(.width == 0.8) %>% 
                    filter(!is.na(time)) %>% 
                    filter(str_ends(name, "_t")),
                  mapping = aes(x = time, y = value, ymin = .lower, ymax = .upper, fill = distribution, color = distribution),
                  alpha = 0.5,
                  key_glyph = "rect") +
  geom_point(data = true_generated_quantities %>%
               filter(!is.na(time)) %>% 
               filter(str_ends(name, "_t")),
             mapping = aes(x = time, y = true_value)) +
  geom_line(data = true_generated_quantities %>%
              filter(!is.na(time)) %>% 
              filter(str_ends(name, "_t")),
            mapping = aes(x = time, y = true_value)) +
  facet_wrap(. ~ name, scales = "free_y",
             labeller = my_labeller_fn) +
  scale_y_continuous(name = "Value", labels = comma) +
  scale_x_continuous(name = "Time") +
  scale_fill_discrete(name = "Distribution", labels = str_to_title) +
  scale_color_discrete(name = "Distribution", labels = str_to_title) +
  ggtitle("Prior and Posterior Credible Intervals for Time-Varying Parameters",
          subtitle = "One simulated dataset, 80% credible intervals, true values in black") +
  theme(legend.position = c(5/6, 1/4))


save_plot_target_asp(filename = path(defense_figures_dir, "simulated_binned_data_plot", ext = "pdf"),
          plot = simulated_binned_data_plot, ncol = 2, nrow = 2)


save_plot_target_asp(filename = path(defense_figures_dir, "single_generated_quantities_simulation_compartment_plot", ext = "pdf"),
          plot = single_generated_quantities_simulation_compartment_plot,
          ncol = 3, nrow = 2, base_height = 2)

save_plot_target_asp(filename = path(defense_figures_dir, "single_generated_quantities_simulation_scalar_plot", ext = "pdf"),
          plot = single_generated_quantities_simulation_scalar_plot,
          ncol = 3, nrow = 3, base_height = 2)

save_plot_target_asp(filename = path(defense_figures_dir, "single_generated_quantities_simulation_time_varying_plot", ext = "pdf"),
          plot = single_generated_quantities_simulation_time_varying_plot,
          ncol = 3, nrow = 2, base_height = 2)
