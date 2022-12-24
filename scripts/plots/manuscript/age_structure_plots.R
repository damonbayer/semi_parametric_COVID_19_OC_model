library(tidyverse)
library(tidybayes)
library(scales)
source("src/plot_functions.R")

latent_curves <-
  read_csv("illustrative_examples/age_structure/data/data.csv") %>%
  select(t, ends_with("_y"), ends_with("_o")) %>%
  pivot_longer(-t,
    names_to = c("name", "population"),
    names_sep = "_(?=(o|y)$)"
  ) %>%
  pivot_wider(names_from = population, values_from = value) %>%
  mutate(total = y + o) %>%
  pivot_longer(cols = c(-t, -name), names_to = "population") %>%
  mutate(population = population %>%
    fct_relevel("y", "o", "total") %>%
    fct_recode(
      General = "y",
      Vulnerable = "o",
      Combined = "total"
    )) %>%
  mutate(type = "latent") %>% 
  drop_na()

true_ifr_dat <-
  latent_curves %>%
  filter(
    name == "I",
    population != "Combined"
  ) %>%
  pivot_wider(names_from = population, values_from = value) %>%
  group_by(t) %>%
  summarize(value = (General * 0.01 + Vulnerable * 0.1) / (General + Vulnerable)) %>%
  rename(time = t) %>% 
  drop_na()

dat <-
  read_csv("illustrative_examples/age_structure/data/data.csv") %>%
  select(t, starts_with("data")) %>%
  pivot_longer(-t) %>%
  mutate(name = str_remove(name, "data_")) %>%
  mutate(type = "observed") %>% 
  drop_na()

latent_curves %>%
  unite(col = grp, name, population, remove = FALSE) %>%
  ggplot(aes(t, value, color = population, group = grp)) +
  facet_wrap(. ~ name, scale = "free_y") +
  geom_line() +
  cowplot::theme_cowplot()

data_ifr_age_structure_plot <-
  ggplot(mapping = aes(t, value)) +
  facet_wrap(. ~ name, scale = "free_y", labeller = as_labeller(. %>% str_replace("_", " ") %>% str_to_title()), ncol = 1) +
  geom_line(
    mapping = aes(color = population),
    data = latent_curves %>%
      filter(str_detect(name, "latent")) %>%
      mutate(name = str_remove(name, "latent_"))
  ) +
  geom_point(data = dat) +
  scale_y_continuous(name = "Count", labels = comma) +
  scale_x_continuous(name = "Time") +
  labs(color = "Population Group") +
  theme_minimal_grid() +
  ggtitle("Latent and Observed Cases and Deaths") +
  theme(legend.position = "bottom")

generated_quantities <-
  read_csv("illustrative_examples/age_structure/data/generated_quantities.csv") %>%
  pivot_longer(-c(iteration, chain), names_to = "name_raw") %>% 
  mutate(name = 
           if_else(
             str_detect(name_raw, "\\[\\d+\\]"), 
             str_extract(name_raw, ".+(?=\\[)"),
             name_raw),
         time = name_raw %>% 
           str_extract("(?<=\\[)\\d+(?=\\])") %>% 
           as.integer()) %>% 
  select(name, time, value) %>% 
  group_by(name, time) %>% 
  median_qi(.width = c(0.5, 0.8, 0.95))

posterior_ifr_age_structure_plot <-
  generated_quantities %>%
  filter(name == "IFR_t") %>%
  ggplot(aes(time, value)) +
  geom_lineribbon(
    mapping = aes(ymin = .lower, ymax = .upper),
    color = brewer_line_color, step = "hv", key_glyph = "rect"
  ) +
  geom_step(data = true_ifr_dat %>%
    filter(time <= generated_quantities %>%
      filter(name == "IFR_t") %>%
      pull(time) %>%
      max()),
    size = 1.5,
    linetype = "dotted") +
  scale_x_continuous("Time") +
  scale_y_continuous("IFR", labels = percent) +
  ggtitle(
    "Posterior Infection Fatality Ratio for Heterogeneous Population",
    "Modelled as Homogeneous Population"
  ) +
  my_theme

save_plot(
  filename = path(figures_dir, "data_ifr_age_structure_plot", ext = "pdf"),
  plot = data_ifr_age_structure_plot,
  ncol = 1,
  nrow = 2
)

save_plot(
  filename = path(figures_dir, "posterior_ifr_age_structure_plot", ext = "pdf"),
  plot = posterior_ifr_age_structure_plot,
  base_asp = 2
)
