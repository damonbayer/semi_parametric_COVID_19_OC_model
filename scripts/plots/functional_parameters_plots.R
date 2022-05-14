library(tidyverse)
source("src/plot_functions.R")

functional_params_dat <- 
  read_csv("data/simulated_data/true_generated_quantities_constant_IFR=false_constant_R0=false_constant_alpha=false_double_IFR_0=false_half_R0_0=false_half_S_0=false_half_alpha_0=false_max_t=42.0_seed=1_use_seroprev=true_use_tests=true.csv") %>% 
  select(-iteration, -chain) %>% 
  pivot_longer(everything(), values_to = "true_value") %>% 
  mutate(time = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  mutate(name = if_else(str_detect(name,"^.+\\["), str_extract(name, "^.+(?=\\[)"), name)) %>% 
  filter(name %in% c("S", "E", "I", "R", "D", "C"))

incidence_dat <- 
  functional_params_dat %>% 
  filter(name %in% c("I", "C")) %>% 
  pivot_wider(names_from = name, values_from = true_value) %>% 
  mutate(delta_C = C - lag(C))

tmp_xlim <- c(5, 35)
tmp_ylim <- functional_params_dat %>% 
  filter(time <=  tmp_xlim[2],
         time >= tmp_xlim[1]) %>% 
  filter(name %in% c("E", "I")) %>% 
  pull(true_value) %>% 
  range()

tmp_asp <- as.numeric((tmp_xlim[2] - tmp_xlim[1]) / (tmp_ylim[2] - tmp_ylim[1]))

functional_parameters_abstract_plot <- 
  functional_params_dat %>% 
  filter(name %in% c("E", "I")) %>% 
  ggplot(aes(time, true_value, color = name)) +
  geom_line(size = 2) +
  my_theme +
  theme_nothing() +
  coord_fixed(ratio = tmp_asp, xlim = tmp_xlim, ylim = tmp_ylim)


N_EI_plot <- 
  ggplot(incidence_dat, aes(time, C)) +
  geom_line() +
  scale_y_continuous(name = TeX("$N_{EI}(t)$"),
                     labels = comma) +
  scale_x_continuous(name = "t") +
  ggtitle(TeX("$N_{EI}(t)$"))

I_plot <- 
  ggplot(incidence_dat, aes(time, I)) +
  geom_line() +
  scale_y_continuous(name = TeX("$I(t)$"),
                     labels = comma) +
  scale_x_continuous(name = "t") +
  ggtitle(TeX("$I(t)$"))


delta_N_EI_plot <- 
  ggplot(incidence_dat %>% drop_na(), aes(time, delta_C)) +
  geom_point() +
  scale_y_continuous(name = TeX("$\\Delta N_{EI}(t_l)$"),
                     labels = comma) +
  scale_x_continuous(name = TeX("$t_l$")) +
  ggtitle(TeX("$\\Delta N_{EI}(t_l)$"))

functional_parameters_ex_plot <- cowplot::plot_grid(I_plot, N_EI_plot, delta_N_EI_plot, nrow = 1, align = "hv")


save_plot_target_asp(filename = "figures/advancement_slides/functional_parameters_ex_plot.pdf",
                     plot = functional_parameters_ex_plot,
                     ncol = 3,
                     base_asp = 1730/650*2, base_height = 3)


save_plot(filename = "figures/advancement_slides/functional_parameters_abstract_plot.pdf",
          plot = functional_parameters_abstract_plot,
          base_asp = 1)
