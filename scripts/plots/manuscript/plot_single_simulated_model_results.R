library(tidyverse)
library(tidybayes)
library(fs)
source("src/plot_functions.R")

real_dat <- read_csv("data/oc_data.csv") %>% 
  filter(time <= 42)

simulated_dat <- read_csv("data/simulated_data/simulated_data_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  slice(1) %>% 
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

true_parameters <-
  read_csv("data/simulated_data/true_parameters_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  select(-iteration, -chain) %>% 
  pivot_longer(everything(), values_to = "true_value") %>% 
  mutate(time = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  mutate(name = if_else(str_detect(name,"^.+\\["), str_extract(name, "^.+(?=\\[)"), name))

true_gq <-
  read_csv("data/simulated_data/true_generated_quantities_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  select(-iteration, -chain) %>% 
  pivot_longer(everything(), values_to = "true_value") %>% 
  mutate(time = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  mutate(name = if_else(str_detect(name,"^.+\\["), str_extract(name, "^.+(?=\\[)"), name))

prior_parameters_summary <- 
  read_csv("results/simulated_posterior_samples_summary/simulated_prior_samples_summary_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  rename(name = parameters,
         .lower = `10.0%`,
         .upper = `90.0%`)

prior_gq_summary <- 
  read_csv("results/simulated_generated_quantities_summary/simulated_prior_generated_quantities_summary_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  rename(name = parameters,
         .lower = `10.0%`,
         .upper = `90.0%`)

single_sim_parameters_summary <- read_csv("results/simulated_posterior_samples_summary/simulated_posterior_samples_summary_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  rename(name = parameters,
         .lower = `10.0%`,
         .upper = `90.0%`) %>% 
  mutate(.width = 0.8)

single_sim_gq_summary <- read_csv("results/simulated_generated_quantities_summary/simulated_generated_quantities_summary_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  rename(name = parameters,
         .lower = `10.0%`,
         .upper = `90.0%`) %>% 
  mutate(.width = 0.8)

all_parameter_intervals <- 
  bind_rows(
    mutate(prior_parameters_summary, dist = "Prior"),
    mutate(single_sim_parameters_summary, dist = "Posterior")) %>% 
  mutate(.width = 0.8) %>%
  mutate(dist = fct_inorder(dist)) %>% 
  select(name, mean, .lower, .upper, dist) %>% 
  mutate(time = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  mutate(name = if_else(str_detect(name,"^.+\\["), str_extract(name, "^.+(?=\\[)"), name))

all_gq_intervals <- 
  bind_rows(
    mutate(prior_gq_summary, dist = "Prior"),
    mutate(single_sim_gq_summary, dist = "Posterior")) %>% 
  mutate(.width = 0.8) %>% 
  mutate(dist = fct_inorder(dist)) %>% 
  select(name, mean, .lower, .upper, dist) %>% 
  mutate(time = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  mutate(name = if_else(str_detect(name,"^.+\\["), str_extract(name, "^.+(?=\\[)"), name))

single_gq_simulation_scalar_plot <- 
  ggplot() +
  geom_pointinterval(data = all_gq_intervals %>%
                       filter(is.na(time)) %>% 
                       filter(!(name %in%  c("S_SEI", "I_EI"))),
                     mapping = aes(y = dist, x = mean, xmin = .lower, xmax = .upper, color = dist),
                     size = line_size) +
  geom_vline(data = true_gq %>%
               filter(is.na(time)) %>% 
               filter(!(name %in%  c("S_SEI", "I_EI"))),
             mapping = aes(xintercept = true_value),
             size = line_size) +
  facet_wrap(. ~ name, scales = "free_x",
             labeller = my_labeller_fn) +
  scale_x_continuous(name = "Value") +
  scale_y_discrete(name = "Distribution") +
  scale_fill_discrete(name = "Distribution") +
  scale_color_discrete(name = "Distribution") +
  ggtitle("Prior and Posterior Credible Intervals for Time-Stationary Parameters",
          subtitle = "One simulated dataset, 80% credible intervals, true values in black") +
  theme(legend.position = "none")


single_gq_simulation_compartment_plot <- 
  ggplot() +
  geom_lineribbon(data = all_gq_intervals %>%
                    filter(!is.na(time)) %>% 
                    filter(name %in% c("S", "E", "I", "R", "D")) %>% 
                    mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))),
                  mapping = aes(x= time, y = mean, ymin = .lower, ymax = .upper, fill = dist, color = dist),
                  alpha = 0.5,
                  key_glyph = "rect") +
  geom_point(data = true_gq %>%
               filter(!is.na(time)) %>% 
               filter(name %in% c("S", "E", "I", "R", "D")) %>% 
               mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))),
             mapping = aes(x = time, y = true_value)) +
  geom_line(data = true_gq %>%
              filter(!is.na(time)) %>% 
              filter(name %in% c("S", "E", "I", "R", "D")) %>% 
              mutate(name = fct_relevel(name, c("S", "E", "I", "R", "D"))),
            mapping = aes(x = time, y = true_value)) +
  facet_wrap(. ~ name, scales = "free_y",
             labeller = my_labeller_fn) +
  scale_y_continuous(name = "Count", labels = comma) +
  scale_x_continuous(name = "Time") +
  scale_fill_discrete(name = "Distribution") +
  scale_color_discrete(name = "Distribution") +
  ggtitle("Prior and Posterior Credible Intervals for Compartments",
          subtitle = "One simulated dataset, 80% credible intervals, true values in black") +
  theme(legend.position = c(5/6, 1/4))

single_gq_simulation_time_varying_plot <- 
  ggplot() +
  geom_lineribbon(data = all_gq_intervals %>%
                    filter(!is.na(time)) %>% 
                    filter(str_ends(name, "_t")),
                  mapping = aes(x= time, y = mean, ymin = .lower, ymax = .upper, fill = dist, color = dist),
                  alpha = 0.5,
                  key_glyph = "rect") +
  geom_point(data = true_gq %>%
               filter(!is.na(time)) %>% 
               filter(str_ends(name, "_t")),
             mapping = aes(x = time, y = true_value)) +
  geom_line(data = true_gq %>%
              filter(!is.na(time)) %>% 
              filter(str_ends(name, "_t")),
            mapping = aes(x = time, y = true_value)) +
  facet_wrap(. ~ name, scales = "free_y",
             labeller = my_labeller_fn) +
  scale_y_continuous(name = "Value", labels = comma) +
  scale_x_continuous(name = "Time") +
  scale_fill_discrete(name = "Distribution") +
  scale_color_discrete(name = "Distribution") +
  ggtitle("Prior and Posterior Credible Intervals for Time-Varying Parameters",
          subtitle = "One simulated dataset, 80% credible intervals, true values in black") +
  theme(legend.position = c(5/6, 1/4))


save_plot(filename = path(figures_dir, "simulated_binned_data_plot", ext = "pdf"),
          plot = simulated_binned_data_plot, ncol = 2, nrow = 2)

save_plot(filename = path(figures_dir, "single_gq_simulation_compartment_plot", ext = "pdf"),
          plot = single_gq_simulation_compartment_plot,
          ncol = 3, nrow = 2)

save_plot(filename = path(figures_dir, "single_gq_simulation_scalar_plot", ext = "pdf"),
          plot = single_gq_simulation_scalar_plot,
          ncol = 3, nrow = 3)

save_plot(filename = path(figures_dir, "single_gq_simulation_time_varying_plot", ext = "pdf"),
          plot = single_gq_simulation_time_varying_plot,
          ncol = 3, nrow = 2)
