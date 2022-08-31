library(tidyverse)
library(fs)
library(tidybayes)

target_model_design <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 1, as.integer(commandArgs(trailingOnly=T[1])))

file_path <- tibble(file_path = dir_ls("results/posterior_predictive//")) %>% 
  mutate(model_design = file_path %>% 
           str_extract("(?<=model_design=)\\d+") %>% 
           as.integer()) %>% 
  filter(model_design == target_model_design) %>% 
  pull(file_path)

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

tidy_predictive_file <- function(file_name) {
  read_csv(file_name) %>% 
    mutate(draw = tidybayes:::draw_from_chain_and_iteration_(chain = chain, iteration = iteration), .after = iteration) %>% 
    pivot_longer(-c(iteration, chain, draw), names_to = "name_raw") %>% 
    mutate(name_raw = str_remove(name_raw, "data_new_|data_")) %>% 
    mutate(name = name_raw %>% str_extract("^.+(?=\\[\\d+\\])"),
           index = name_raw%>% str_extract("(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
    select(-name_raw) %>% 
    pivot_wider(values_from = value) %>% 
    left_join(dat %>% select(index, tests)) %>% 
    mutate(test_positivity = cases / tests) %>%
    select(-tests) %>% 
    pivot_longer(-c(iteration, chain, draw, index)) %>% 
    rename_with(~str_c(".", .), c(iteration, chain, draw)) %>% 
    drop_na()
}

prep_predictive_for_plotting <- function(tidy_predictive, index_date_conversion) {
  tidy_predictive %>% 
    select(index, name, value) %>% 
    group_by(index, name) %>% 
    median_qi(.width = ci_widths) %>% 
    left_join(index_date_conversion) %>% 
    mutate(date = if_else(str_detect(name, "seroprev"), oc_seroprev_data$end_date, date)) %>% 
    select(date, everything(), -index) %>% 
    mutate(name = fct_inorder(name))
}

score_predictions <- function(tidy_predictive, index_date_conversion) {
  tidy_predictive %>% 
    filter(name != "test_positivity") %>% 
    group_by(index, name) %>% 
    summarize(pp_draws = list(value),
              .groups = "drop") %>% 
    left_join(dat %>%
                select(index, date, cases, deaths) %>% 
                pivot_longer(-c(index, date), values_to = "obs_value")) %>% 
    mutate(ll = map2_dbl(obs_value, pp_draws, ~log(mean(.x == .y)))) %>% 
    left_join(index_date_conversion) %>% 
    select(date, name, ll) %>% 
    pivot_wider(names_from = name, values_from = ll, names_prefix = "ll_")
}

tidy_predictive <- tidy_predictive_file(file_path)

pp_for_plotting <- prep_predictive_for_plotting(tidy_predictive = tidy_predictive,
                                                index_date_conversion = index_date_conversion)

prediction_score <- score_predictions(tidy_predictive = tidy_predictive,
                                      index_date_conversion = index_date_conversion)

dir_create("results/tidy_posterior_predictive")
write_csv(pp_for_plotting, str_replace_all(file_path, "posterior_predictive", "tidy_posterior_predictive"))

dir_create("results/prediction_score")
write_csv(prediction_score, str_replace_all(file_path, "posterior_predictive", "prediction_score"))