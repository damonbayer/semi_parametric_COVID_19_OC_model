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
  mutate(across(c(wall, compute), ~seconds(.) |> seconds_to_period()))


duration_dat |> 
  left_join(model_table) |> 
  arrange(desc(compute))

duration_dat |> 
  left_join(model_table) |> 
  filter(constant_alpha == F,
         constant_IFR == F,
         constant_R0 == F,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F,
         max_t == 42,
         use_seroprev == T,
         use_tests == T) |> 
  select(wall, compute)
  
