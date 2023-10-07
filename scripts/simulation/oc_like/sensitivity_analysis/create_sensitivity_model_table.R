# Setup design for all models fit in paper
library(tidyverse)
library(fs)
model_table <-
  diag(4) %>% 
  as_tibble() %>%
  mutate(across(everything(), ~!as.logical(.))) %>% 
  `colnames<-`(c("use_deaths", "use_tests", "use_seroprev", "use_wide_priors")) %>% 
  mutate(use_wide_priors = !use_wide_priors) %>% 
  add_row(.before = 1, tibble(use_deaths = TRUE, use_tests = TRUE, use_seroprev = TRUE, use_wide_priors = FALSE)) %>% 
  mutate(model_design = row_number()) %>% 
  expand_grid(data_id= 1:200) %>%
  mutate(sim_id = row_number()) %>% 
  select(sim_id, model_design, data_id, everything())

write_csv(model_table, path("scripts", "simulation", "oc_like", "sensitivity_analysis", "sensitivity_model_table", ext = "csv"))



model_table %>% 
  filter(model_design >= 2) %>% 
  summarize(min(sim_id), max(sim_id))
unique(model_table$model_design)
