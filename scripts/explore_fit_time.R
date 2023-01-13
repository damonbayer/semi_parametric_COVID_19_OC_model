library(tidyverse)
library(fs)
library(lubridate)

model_table <- read_csv("model_table.csv") %>% distinct(model_design, .keep_all = T)

duration_dat <- 
  imap_dfr(dir_ls("results/duration"), ~mutate(read_csv(.x), full_path = .y)) |> 
  relocate(full_path) |> 
  mutate(model_id = full_path %>% 
           str_extract("(?<=model_id=)\\d+") %>% 
           as.numeric()) |> 
  mutate(across(c(wall, compute), ~seconds(.) |> seconds_to_period())) %>% 
  select(model_id, wall, compute) %>% 
  left_join(model_table)


oc_lik_sim_duration_dat <- 
  imap_dfr(dir_ls("results/simulation/oc_like/duration"), ~mutate(read_csv(.x), full_path = .y)) %>% 
  mutate(sim_id = full_path %>% path_ext_remove() %>% str_extract("\\d+$") %>% as.integer()) %>% 
  arrange(sim_id) %>% 
  mutate(across(c(wall, compute), ~seconds(.) |> seconds_to_period())) %>% 
  select(sim_id, wall, compute)

oc_lik_sim_duration_dat %>% 
  mutate(across(c(wall, compute), period_to_seconds)) %>% 
  summarize(across(c(wall, compute), mean)) %>% 
  mutate(across(c(wall, compute), ~seconds(.) |> seconds_to_period()))