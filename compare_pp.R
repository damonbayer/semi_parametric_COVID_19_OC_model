library(tidyverse)
library(tidybayes)
library(cowplot)
library(scales)
ci_widths <- c(0.5, 0.8, 0.95)
# NOTE THIS FILE DOES NOT HANDLE SEROPREV TIMES CORRECTLY
# SHOULD SAVE COMMON AESTHETICS IN A NEW THEME
dat <- 
  read_csv("data/oc_data.csv") %>% 
  mutate(index = 1:n()) %>% 
  rename(date = end_date) %>% 
  select(-start_date)

dat_tidy <-
  dat %>% 
  select(date, cases, tests, deaths) %>% 
  mutate(test_positivity = cases / tests) %>% 
  select(-tests) %>% 
  pivot_longer(-date)

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
    mutate(name = fct_inorder(name))
  
  vector_intervals <- 
    tidy_gq %>% 
    filter(!is.na(index)) %>% 
    select(-starts_with(".")) %>% 
    group_by(name, index) %>% 
    median_qi(.width = ci_widths) %>% 
    left_join(index_date_conversion) %>% 
    select(date, everything(), -index) %>% 
    mutate(name = fct_inorder(name))
  
  list(scalar_samples = scalar_samples,
       vector_intervals = vector_intervals)
}

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
    select(date, everything(), -index) %>% 
    mutate(name = fct_inorder(name)) 
}

gq_forecast_randn <-
  tidy_gq_file("gq_forecast_randn.csv") %>% 
  prep_gq_for_plotting(index_date_conversion)

gq_forecast_zeros <-
  tidy_gq_file("gq_forecast_zeros.csv") %>% 
  prep_gq_for_plotting(index_date_conversion)


predictive_randn <- 
  tidy_predictive_file("predictive_randn.csv") %>% 
  prep_predictive_for_plotting(index_date_conversion)

predictive_zeros <- 
  tidy_predictive_file("predictive_zeros.csv") %>% 
  prep_predictive_for_plotting(index_date_conversion)


gq_forecast_zeros$vector_intervals %>% 
  ggplot(aes(date, value, ymin = .lower, ymax = .upper)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon() +
  scale_y_continuous(label = comma) +
  scale_fill_brewer(label = ~percent(as.numeric(.)), guide = guide_legend(reverse = TRUE)) +
  theme_minimal_grid()

gq_forecast_zeros$scalar_samples %>% 
  ggplot(aes(value)) +
  facet_wrap(. ~ name, scales = "free_x") +
  stat_halfeye(normalize = "panels") +
  theme_minimal()


predictive_randn %>% 
  ggplot(aes(date, value)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon(mapping = aes(ymin = .lower, ymax = .upper)) +
  geom_point(data = dat_tidy %>% filter(date <= max(predictive_randn$date))) +
  scale_y_continuous(label = comma) +
  scale_fill_brewer(label = ~percent(as.numeric(.)), guide = guide_legend(reverse = TRUE)) +
  theme_minimal_grid()

predictive_zeros %>% 
  ggplot(aes(date, value)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_lineribbon(mapping = aes(ymin = .lower, ymax = .upper)) +
  geom_point(data = dat_tidy %>% filter(date <= max(predictive_zeros$date))) +
  scale_y_continuous(label = comma) +
  scale_fill_brewer(label = ~percent(as.numeric(.)), guide = guide_legend(reverse = TRUE)) +
  theme_minimal_grid()

