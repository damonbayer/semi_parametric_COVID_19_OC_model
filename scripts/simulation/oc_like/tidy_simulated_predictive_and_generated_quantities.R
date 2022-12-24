library(tidyverse)
library(fs)
library(tidybayes)
library(posterior)
library(scoringutils)

sim_id <- ifelse(length(commandArgs(trailingOnly=T)) == 0, 1, as.integer(commandArgs(trailingOnly=T[1])))
sim_name <- "oc_like"

sim_results_path <- path("results", "simulation", sim_name)

target_max_t <- 42

ci_widths <- c(0.5, 0.8, 0.95)

create_draws <- function(file_path) {
  read_csv(file_path) %>%
    rename_with(~str_c(".", .), c(iteration, chain)) %>%
    as_draws()
}

tidy_generated_quantities <- function(generated_quantities_draws) {
  generated_quantities_draws %>%
    select(-any_of("lp")) %>% 
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
      name == "seroprev_mean" ~ 20L,
      TRUE ~ index - 1L)) %>% 
    select(name, time, value, starts_with(".")) %>% 
    arrange(name, time, .width)
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
    mutate(time = if_else(name == "seroprev_cases", 20L, index)) %>% 
    select(name, time, value, starts_with(".")) %>% 
    arrange(name, time, .width)
}

compute_generated_quantities_sd <- function(generated_quantities_draws) {
  generated_quantities_draws %>%
    select(-any_of("lp")) %>% 
    summarize_draws(sd,
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
      name == "seroprev_mean" ~ 20L,
      TRUE ~ index - 1L)) %>% 
    select(name, time, sd) %>% 
    arrange(name, time)
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
      name == "seroprev_mean" ~ 20L,
      TRUE ~ index - 1L)) %>% 
    select(name, time, everything(), -index, -time)
}

if (sim_id == 0) {
  prior_generated_quantities_file_path <- path(sim_results_path, "prior_generated_quantities", "prior_generated_quantities", ext = "csv")
  prior_predictive_file_path <- path(sim_results_path, "prior_predictive", "prior_predictive", ext = "csv")
  
  prior_generated_quantities_draws <- create_draws(prior_generated_quantities_file_path)
  tidy_prior_generated_quantities <- tidy_generated_quantities(prior_generated_quantities_draws)
  prior_generated_quantities_sd <- compute_generated_quantities_sd(prior_generated_quantities_draws)

  prior_predictive_draws <- create_draws(prior_predictive_file_path)
  tidy_prior_predictive <- tidy_predictive(prior_predictive_draws)
  
  dir_create(path(sim_results_path, "tidy_prior_generated_quantities"))
  dir_create(path(sim_results_path, "prior_generated_quantities_sd"))
  dir_create(path(sim_results_path, "tidy_prior_predictive"))

  write_csv(prior_generated_quantities_sd, path(sim_results_path, "prior_generated_quantities_sd", "prior_generated_quantities_sd", ext = "csv"))
  write_csv(tidy_prior_generated_quantities, path(sim_results_path, "tidy_prior_generated_quantities", "tidy_prior_generated_quantities", ext = "csv"))
  write_csv(tidy_prior_predictive, path(sim_results_path, "tidy_prior_predictive", "tidy_prior_predictive", ext = "csv"))
} else {
  posterior_generated_quantities_file_path <- path("results/simulation/oc_like/posterior_generated_quantities", str_c("posterior_generated_quantities_sim_id=", sim_id), ext = "csv")
  posterior_predictive_file_path <- path("results/simulation/oc_like/posterior_predictive", str_c("posterior_predictive_sim_id=", sim_id), ext = "csv")
  
  posterior_generated_quantities_draws <- create_draws(posterior_generated_quantities_file_path)
  tidy_posterior_generated_quantities <- tidy_generated_quantities(posterior_generated_quantities_draws)
  posterior_generated_quantities_sd <- compute_generated_quantities_sd(posterior_generated_quantities_draws)
  
  posterior_lp <- extract_lp_tbl(posterior_generated_quantities_draws)
  posterior_diagnostics <- compute_diagnostics(posterior_generated_quantities_draws)
  
  posterior_predictive_draws <- create_draws(posterior_predictive_file_path)
  tidy_posterior_predictive <- tidy_predictive(posterior_predictive_draws)
  
  dir_create(path(sim_results_path, "tidy_posterior_generated_quantities"))
  dir_create(path(sim_results_path, "posterior_generated_quantities_sd"))
  dir_create(path(sim_results_path, "posterior_lp"))
  dir_create(path(sim_results_path, "posterior_diagnostics"))
  dir_create(path(sim_results_path, "tidy_posterior_predictive"))
  
  write_csv(tidy_posterior_generated_quantities, path(sim_results_path, "tidy_posterior_generated_quantities", str_c("tidy_posterior_generated_quantities_sim_id=", sim_id), ext = "csv"))
  write_csv(posterior_generated_quantities_sd, path(sim_results_path, "posterior_generated_quantities_sd", str_c("posterior_generated_quantities_sd_sim_id=", sim_id), ext = "csv"))
  write_csv(posterior_lp, path(sim_results_path, "posterior_lp", str_c("posterior_lp_sim_id=", sim_id), ext = "csv"))
  write_csv(posterior_diagnostics, path(sim_results_path, "posterior_diagnostics", str_c("posterior_diagnostics_sim_id=", sim_id), ext = "csv"))
  write_csv(tidy_posterior_predictive, path(sim_results_path, "tidy_posterior_predictive", str_c("tidy_posterior_predictive_sim_id=", sim_id), ext = "csv"))
}
