library(tidyverse)
library(tidybayes)
library(scales)
source("src/plot_functions.R")

latent_curves <- 
  read_csv("illustrative_examples/age_structure/data/data.csv") %>% 
  select(t, ends_with("_y"), ends_with("_o")) %>% 
  pivot_longer(-t,
               names_to = c("name", "age_group"),
               names_sep = "_(?=(o|y)$)") %>% 
  pivot_wider(names_from = age_group, values_from = value) %>% 
  mutate(total = y + o) %>% 
  pivot_longer(cols = c(-t, - name), names_to = "age_group") %>% 
  mutate(age_group = age_group %>% 
           fct_relevel("y", "o", "total") %>% 
           fct_recode(Young = "y",
                      Old = "o",
                      Combined = "total")) %>% 
  mutate(type = "latent")

latent_curves %>% 
  filter(name == "I",
         age_group != "Combined") %>% 
  pivot_wider(names_from = age_group, values_from = value) %>% 
  group_by(t) %>% 
  summarize(true_ifr = (Young * 0.01 + Old * 0.1) / (Young + Old)) %>% 
  ggplot(aes(t, true_ifr)) +
  geom_line()

unique(latent_curves$name)

dat <- 
  read_csv("illustrative_examples/age_structure/data/data.csv") %>% 
  select(t, starts_with("data")) %>% 
  pivot_longer(-t) %>% 
  mutate(name = str_remove(name, "data_")) %>% 
  mutate(type = "observed")

latent_curves %>% 
  unite(col = grp, name, age_group, remove = F) %>% 
  ggplot(aes(t, value, color = age_group, group = grp)) +
  facet_wrap(. ~ name, scale = "free_y") +
  geom_line() +
  cowplot::theme_cowplot()
  
data_ifr_age_structure_plot <-     
  ggplot(mapping = aes(t, value)) +
  facet_wrap(. ~ name, scale = "free_y", labeller = as_labeller(. %>% str_replace("_", " ") %>% str_to_title())) +
  geom_line(mapping = aes(color = age_group),
            data = latent_curves %>% 
              filter(str_detect(name, "latent")) %>% 
              mutate(name = str_remove(name, "latent_"))) +
  geom_point(data = dat) +
  scale_y_continuous(name = "Count", labels = comma) +
  scale_x_continuous(name = "Time") +
  labs(color = "Age Group") +
  theme_minimal_grid() +
  ggtitle("Latent and Observed Cases and Deaths") +
  theme(legend.position = "bottom")

posterior_predictive <- 
  read_csv("illustrative_examples/age_structure/data/posterior_predictive.csv") %>% 
  mutate(draw = tidybayes:::draw_from_chain_and_iteration_(chain = chain, iteration = iteration), .after = iteration) %>% 
  pivot_longer(-c(iteration, chain, draw), names_to = "name_raw") %>% 
  mutate(name_raw = str_remove(name_raw, "data_new_|data_")) %>% 
  mutate(name = name_raw %>% str_extract("^.+(?=\\[\\d+\\])"),
         index = name_raw%>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  relocate(chain, iteration, draw, name, index, value) %>% 
  rename_with(~str_c(".", .), c(iteration, chain, draw))

generated_quantities <-
  read_csv("illustrative_examples/age_structure/data/generated_quantities.csv")  %>% 
  mutate(draw = tidybayes:::draw_from_chain_and_iteration_(chain = chain, iteration = iteration), .after = iteration) %>% 
  pivot_longer(-c(iteration, chain, draw), names_to = "name_raw") %>% 
  separate(col = name_raw,
           into = c("name", "index"),
           sep = "\\[|\\]",
           remove = T,
           fill = "right",
           extra = "drop",
           convert = T) %>% 
  select(chain, iteration, everything()) %>% 
  rename_with(~str_c(".", .), c(iteration, chain, draw)) %>% 
  select(-starts_with(".")) %>% 
  group_by(name, index) %>% 
  median_qi(.width = c(0.5, 0.8, 0.95))


posterior_ifr_age_structure_plot <- 
  generated_quantities %>% 
  filter(name == "IFR_t") %>% 
  ggplot(aes(index, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  scale_x_continuous("Time") +
  scale_y_continuous("IFR", labels = percent, breaks = c(0, 0.01, 0.05, 0.1, 0.15)) +
  annotate("text", x = 2, y = 0.1, label = "True IFR for Old Population", vjust = -0.5, hjust = 0) +
  annotate("text", x = 0, y = 0.01, label = "True IFR for Young Population", vjust = -0.5, hjust = 0) +
  geom_hline(yintercept = c(0.01, 0.1), linetype = "dashed") +
  ggtitle("Posterior Infection Fatality Ratio for Heterogeneous Population",
          "Modelled as Homogeneous Population") +
  my_theme

save_plot_target_asp(filename = "figures/advancement_slides/data_ifr_age_structure_plot.pdf",
                     plot = data_ifr_age_structure_plot,
                     ncol = 2,
                     nrow = 1,
                     base_asp = 16/9)

save_plot_target_asp(filename = "figures/advancement_slides/posterior_ifr_age_structure_plot.pdf",
                     plot = posterior_ifr_age_structure_plot,
                     ncol = 1,
                     nrow = 1,
                     base_height = 5,
                     base_asp = 16/9)
