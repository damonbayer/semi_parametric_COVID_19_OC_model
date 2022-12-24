library(tidyverse)
library(fs)
library(tidybayes)
library(posterior)
library(scoringutils)

target_model_design <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 79, as.integer(commandArgs(trailingOnly=T[1])))

ci_widths <- c(0.5, 0.8, 0.95)

dat <-
  read_csv("data/oc_data.csv") %>%
  select(time, cases, deaths, tests, date = end_date) %>%
  mutate(test_positivity = cases / tests)

dat_tidy <- 
  dat %>% 
  select(-date, -tests, -test_positivity) %>% 
  pivot_longer(-time)

time_interval_in_days <- as.numeric(dat$date[2] - dat$date[1])

time_date_key <- 
  dat %>% 
  select(time, date) %>% 
  add_row(., time = 0, date = min(.$date) - time_interval_in_days, .before = 1) 

oc_seroprev_data <- read_csv("data/oc_seroprev_data.csv")

create_draws <- function(file_path) {
  read_csv(file_path) %>%
    rename_with(~str_c(".", .), c(iteration, chain)) %>%
    as_draws()
}

tidy_generated_quantities <- function(generated_quantities_draws) {
  generated_quantities_draws %>%
    select(-lp) %>% 
    pivot_longer(-starts_with(".")) %>%
    group_by(name) %>%
    median_qi(.width = ci_widths) %>%
    separate(col = name,
             into = c("name", "index"),
             sep = "\\[|\\]",
             remove = T,
             fill = "right",
             extra = "drop",
             convert = T) %>%
    mutate(time = case_when(
      name %in% c("α_t", "cases_bb_mean", "cases_nb_mean", "cases_mean", "deaths_mean") ~ index,
      name == "seroprev_mean" ~ as.integer(oc_seroprev_data$time),
      TRUE ~ index - 1L)) %>% 
    left_join(time_date_key) %>% 
    select(name, date, value, starts_with(".")) %>% 
    arrange(name, date, .width)
}

tidy_predictive <- function(predictive_draws) {
  predictive_draws %>%
    pivot_longer(-starts_with(".")) %>% 
    group_by(name) %>%
    median_qi(.width = ci_widths) %>%
    separate(col = name,
             into = c("name", "index"),
             sep = "\\[|\\]",
             remove = T,
             fill = "right",
             extra = "drop",
             convert = T) %>% 
    mutate(name = str_remove(name, "data_new_|data_")) %>% 
    mutate(time = if_else(name == "seroprev_cases", as.integer(oc_seroprev_data$time), index)) %>% 
    left_join(time_date_key) %>% 
    select(name, date, value, starts_with(".")) %>% 
    arrange(name, date, .width)
}

extract_lp_tbl <- function(generated_quantities_draws) {
  generated_quantities_draws %>% 
    select(starts_with("."), lp) %>% 
    as_tibble() %>% 
    select(-.draw)
}

compute_diagnostics <- function(generated_quantities_draws) {
  summarize_draws(generated_quantities_draws,
                  rhat, rhat_basic, ess_basic, ess_bulk, ess_tail,
                  .cores = parallelly::availableCores()) %>% 
    separate(col = variable,
             into = c("name", "index"),
             sep = "\\[|\\]",
             remove = T,
             fill = "right",
             extra = "drop",
             convert = T) %>% 
    mutate(time = case_when(
      name %in% c("α_t", "cases_bb_mean", "cases_nb_mean", "cases_mean", "deaths_mean") ~ index,
      name == "seroprev_mean" ~ as.integer(oc_seroprev_data$time),
      TRUE ~ index - 1L)) %>% 
    left_join(time_date_key) %>% 
    select(name, date, everything(), -index, -time)
}

score_with_override <- function (data, metrics = NULL, override = NULL, ...) {
  check_data <- check_forecasts(data)
  data <- check_data$cleaned_data
  prediction_type <- check_data$prediction_type
  forecast_unit <- check_data$forecast_unit
  target_type <- check_data$target_type
  metrics <- scoringutils:::check_metrics(metrics)
  
  if (!is.null(override)) {
    prediction_type <- override
    target_type <- override
  }
  
  if (target_type == "binary") {
    scores <- scoringutils:::score_binary(data = data, forecast_unit = forecast_unit, 
                                          metrics = metrics)
  }
  if (prediction_type == "quantile") {
    scores <- scoringutils:::score_quantile(data = data, forecast_unit = forecast_unit, 
                                            metrics = metrics, ...)
  }
  if (prediction_type %in% c("integer", "continuous") && (target_type != 
                                                          "binary")) {
    scores <- scoringutils:::score_sample(data = data, forecast_unit = forecast_unit, 
                                          metrics = metrics, prediction_type = prediction_type)
  }
  return(scores)
}

score_predictions <- function(predictive_draws) {
  predictive_draws %>% 
    pivot_longer(-starts_with(".")) %>% 
    group_by(name) %>%
    separate(col = name,
             into = c("name", "index"),
             sep = "\\[|\\]",
             remove = T,
             fill = "right",
             extra = "drop",
             convert = T) %>% 
    mutate(name = str_remove(name, "data_new_|data_")) %>% 
    mutate(time = if_else(name == "seroprev_cases", as.integer(oc_seroprev_data$time), index)) %>% 
    select(sample = .draw, prediction = value, name, time) %>% 
    mutate(weeks_ahead = time - target_max_t) %>% 
    filter(weeks_ahead > 0) %>% 
    left_join(dat_tidy %>% 
                rename(true_value = value)) %>% 
    left_join(time_date_key) %>% 
    mutate(model = target_model_design) %>% 
    select(date, target_type = name, model, weeks_ahead, prediction, sample, true_value) %>% 
    score_with_override(override = "continuous")
}

model_info <-
  read_csv("model_table.csv") %>% 
  distinct(model_design) %>% 
  left_join(
    tibble(posterior_generated_quantities_file_path = dir_ls("results/posterior_generated_quantities/")) %>% 
      mutate(model_design = posterior_generated_quantities_file_path %>% 
               str_extract("(?<=model_design=)\\d+") %>% 
               as.integer())
  ) %>%
  left_join(
    tibble(prior_generated_quantities_file_path = dir_ls("results/prior_generated_quantities/")) %>% 
      mutate(model_design = prior_generated_quantities_file_path %>% 
               str_extract("(?<=model_design=)\\d+") %>% 
               as.integer())
  ) %>%
  left_join(
    tibble(posterior_predictive_file_path = dir_ls("results/posterior_predictive/")) %>% 
      mutate(model_design = posterior_predictive_file_path %>% 
               str_extract("(?<=model_design=)\\d+") %>% 
               as.integer())
  ) %>%
  left_join(
    tibble(prior_predictive_file_path = dir_ls("results/prior_predictive/")) %>% 
      mutate(model_design = prior_predictive_file_path %>% 
               str_extract("(?<=model_design=)\\d+") %>% 
               as.integer())
  ) %>% 
  filter(model_design == target_model_design)

prior_generated_quantities_file_path <- model_info$prior_generated_quantities_file_path
posterior_generated_quantities_file_path <- model_info$posterior_generated_quantities_file_path
prior_predictive_file_path <- model_info$prior_predictive_file_path
posterior_predictive_file_path <- model_info$posterior_predictive_file_path

target_max_t <- posterior_predictive_file_path %>% str_extract("(?<=max_t=)\\d+") %>% as.integer()

prior_generated_quantities_draws <- create_draws(prior_generated_quantities_file_path)
posterior_generated_quantities_draws <- create_draws(posterior_generated_quantities_file_path)

tidy_prior_generated_quantities <- tidy_generated_quantities(prior_generated_quantities_draws)
tidy_posterior_generated_quantities <- tidy_generated_quantities(posterior_generated_quantities_draws)

posterior_lp <- extract_lp_tbl(posterior_generated_quantities_draws)
posterior_diagnostics <- compute_diagnostics(posterior_generated_quantities_draws)

prior_predictive_draws <- create_draws(prior_predictive_file_path)
posterior_predictive_draws <- create_draws(posterior_predictive_file_path)

tidy_prior_predictive <- tidy_predictive(prior_predictive_draws)
tidy_posterior_predictive <- tidy_predictive(posterior_predictive_draws)

posterior_predictive_score <- score_predictions(posterior_predictive_draws)

dir_create("results/tidy_prior_generated_quantities")
dir_create("results/tidy_posterior_generated_quantities")
dir_create("results/posterior_lp")
dir_create("results/posterior_diagnostics")
dir_create("results/tidy_prior_predictive")
dir_create("results/tidy_posterior_predictive")
dir_create("results/posterior_predictive_score")

write_csv(tidy_prior_generated_quantities, str_replace_all(prior_generated_quantities_file_path, "prior_generated_quantities", "tidy_prior_generated_quantities"))
write_csv(tidy_posterior_generated_quantities, str_replace_all(posterior_generated_quantities_file_path, "posterior_generated_quantities", "tidy_posterior_generated_quantities"))
write_csv(posterior_lp, str_replace_all(posterior_generated_quantities_file_path, "posterior_generated_quantities", "posterior_lp"))
write_csv(posterior_diagnostics, str_replace_all(posterior_generated_quantities_file_path, "posterior_generated_quantities", "posterior_diagnostics"))
write_csv(tidy_prior_predictive, str_replace_all(prior_predictive_file_path, "prior_predictive", "tidy_prior_predictive"))
write_csv(tidy_posterior_predictive, str_replace_all(posterior_predictive_file_path, "posterior_predictive", "tidy_posterior_predictive"))
write_csv(posterior_predictive_score, str_replace_all(posterior_predictive_file_path, "posterior_predictive", "posterior_predictive_score"))
