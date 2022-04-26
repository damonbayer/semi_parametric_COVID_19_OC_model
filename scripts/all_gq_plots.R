library(gridExtra)
library(tidyverse)
library(tidybayes)
library(fs)
library(scales)
library(latex2exp)
library(cowplot)

model_table <- read_csv("model_table.csv")

vars_to_color <- 
  model_table %>% 
  select(where(~is.logical(.))) %>% 
  colnames()

read_gq_dir <- function(path) {
  tibble(full_path = dir_ls(path)) %>% 
    mutate(data = map(full_path, read_csv)) %>% 
    mutate(model_id = full_path %>% 
             str_extract("(?<=model_id=)\\d+") %>% 
             as.numeric()) %>% 
    left_join(model_table %>% distinct(model_design, .keep_all = T)) %>% 
    select(-full_path, -model_id, -seed) %>% 
    unnest(data) %>% 
    mutate(plot_name = case_when(
      name == "cases_bb_mean" ~ "Cases BB Mean",
      name == "cases_mean" ~ "Cases Mean",
      name == "cases_nb_mean" ~ "Cases NB Mean",
      name == "D" ~ "D",
      name == "deaths_mean" ~ "Deaths Mean",
      name == "dur_infectious_days" ~ "Infectious Period (Days)",
      name == "dur_latent_days" ~ "Latent Period (Days)",
      name == "E" ~ "E",
      name == "I" ~ "I",
      name == "I_EI" ~ "Init I / (E + I)",
      name == "IFR_t" ~ "IFR",
      name == "R" ~ "R",
      name == "R₀_t" ~ "R_0",
      name == "Rₜ_t" ~ "R_t",
      name == "S" ~ "S",
      name == "S_SEI" ~ "Init S / (S + E + I)",
      name == "seroprev_mean" ~ "Seroprev Mean",
      name == "α_t" ~ "\\alpha",
      name == "β_t" ~ "\\beta",
      name == "ρ_cases_t" ~ "\\rho Cases",
      name == "ρ_death" ~ "\\rho Deaths",
      name == "σ_IFR" ~ "\\sigma IFR",
      name == "σ_R0" ~ "\\sigma R_0",
      name == "σ_α" ~ "\\sigma \\alpha",
      name == "σ_ρ_cases" ~ "\\sigma \\rho Cases",
      name == "ϕ_cases_bb" ~ "\\phi Cases BB",
      name == "ϕ_cases_nb" ~ "\\phi Cases NB",
      name == "ϕ_deaths" ~ "\\phi Deaths")) %>% 
    mutate(plot_name = fct_relabel(plot_name, ~(TeX(., output = "character"))))
}

all_vector_gq <- read_gq_dir("results/tidy_vector_generated_quantities/")
all_scalar_gq <- read_gq_dir("results/tidy_scalar_generated_quantities/")
 
all_vector_gq <- bind_rows(
  all_vector_gq,
  all_scalar_gq %>% 
    filter(name %in% c("α", "β", "R₀", "IFR")) %>% 
    mutate(name = str_c(name, "_t")) %>% 
    right_join(distinct(all_vector_gq, model_design, date)) %>% 
    drop_na()
  )

all_scalar_gq <-
  all_scalar_gq %>% 
  filter(!(name %in% c("α", "β", "R₀", "IFR")))

# Time Varying Plots ------------------------------------------------------

make_time_varying_plot <- function(var_to_color) {
  all_vector_gq %>% 
    filter(max_t == 42) %>% 
    filter(.width == 0.95) %>% 
    ggplot(aes_string("date", "value", ymin = ".lower", ymax = ".upper", group = "model_design", color = var_to_color)) +
    facet_wrap(. ~ plot_name, scales = "free_y", labeller = label_parsed) +
    geom_line(size = 1, alpha = 0.5) +
    theme_cowplot() +
    theme(legend.position = "bottom")
}

all_time_varying_plots <- map(vars_to_color, make_time_varying_plot)


 # Scalar Plots ------------------------------------------------------------
make_scalar_plot <- function(var_to_color) {
  all_scalar_gq %>% 
    filter(max_t == 42) %>% 
    mutate(model_design = as_factor(model_design)) %>% 
    ggplot(aes_string("value", "model_design", xmin = ".lower", xmax = ".upper", color = var_to_color)) +
    facet_wrap(. ~ plot_name, scales = "free_x", labeller =  label_parsed) +
    geom_pointinterval() +
    theme_cowplot() +
    theme(legend.position = "bottom")
}

all_scalar_plots <- map(vars_to_color, make_scalar_plot)

# Save Plots --------------------------------------------------------------
dir_create("figures")

ggsave2(filename = "figures/all_time_varying_plots.pdf",
        plot = marrangeGrob(all_time_varying_plots, nrow=1, ncol=1),
        width = 8,
        height = 6,
        scale = 2)

ggsave2(filename = "figures/all_scalar_plots.pdf",
        plot = marrangeGrob(all_scalar_plots, nrow=1, ncol=1),
        width = 8,
        height = 8)
