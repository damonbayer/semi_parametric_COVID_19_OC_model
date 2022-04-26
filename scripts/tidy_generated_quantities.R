library(tidyverse)
library(fs)
library(tidybayes)

model_design <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 1, as.integer(commandArgs(trailingOnly=T[1])))
file_path <- tibble(full_path = dir_ls("results/generated_quantities/")) %>% 
  mutate(obs_model_design = full_path %>% 
           str_extract("(?<=model_design=)\\d+") %>%
           as.integer()) %>% 
  filter(obs_model_design == model_design) %>% 
  slice(1) %>% 
  pull(full_path)

ci_widths <- c(0.5, 0.8, 0.95)

dat <- 
  read_csv("data/oc_data.csv") %>% 
  mutate(index = 1:n()) %>% 
  rename(date = end_date) %>% 
  select(-start_date)

max_data_date <- 
  dat %>% 
  filter(time == 46) %>% 
  pull(date)

dat_tidy <-
  dat %>% 
  select(date, cases, tests, deaths) %>% 
  mutate(test_positivity = cases / tests) %>% 
  select(-tests) %>% 
  pivot_longer(-date)

oc_seroprev_data <- read_csv("data/oc_seroprev_data.csv")

time_interval_in_days <- as.numeric(dat$date[2] - dat$date[1])

index_date_conversion <- 
  dat %>% 
  select(date, index) %>% 
  add_row(date = dat$date[1] - time_interval_in_days, index = 0, .before = 1)

tidy_gq_file <- function(file_name) {
  read_csv(file_name) %>% 
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
    rename_with(~str_c(".", .), c(iteration, chain, draw))
}

prep_gq_for_plotting <- function(tidy_gq, index_date_conversion) {
  scalar_samples <- 
    tidy_gq %>% 
    filter(is.na(index)) %>% 
    select(-index) %>% 
    mutate(name = fct_inorder(name)) %>% 
    select(-starts_with(".")) %>% 
    group_by(name) %>% 
    median_qi(.width = ci_widths)
  
  vector_intervals <- 
    tidy_gq %>% 
    filter(!is.na(index)) %>% 
    select(-starts_with(".")) %>% 
    group_by(name, index) %>% 
    median_qi(.width = ci_widths) %>% 
    left_join(index_date_conversion) %>% 
    select(date, everything(), -index) %>% 
    mutate(name = fct_inorder(name)) %>% 
    mutate(date = if_else(str_detect(name, "seroprev"), oc_seroprev_data$end_date, date))
  
  list(scalar_samples = scalar_samples,
       vector_intervals = vector_intervals)
}

gq_for_plotting <-
  tidy_gq_file(file_path) %>% 
  prep_gq_for_plotting(index_date_conversion = index_date_conversion)


dir_create("results/tidy_scalar_generated_quantities")
dir_create("results/tidy_vector_generated_quantities")

write_csv(gq_for_plotting$scalar_samples, str_replace_all(file_path, "generated_quantities", "tidy_scalar_generated_quantities"))
write_csv(gq_for_plotting$vector_intervals, str_replace_all(file_path, "generated_quantities", "tidy_vector_generated_quantities"))
