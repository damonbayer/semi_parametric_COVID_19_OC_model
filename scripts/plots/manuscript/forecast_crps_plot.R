library(tidyverse)
library(fs)
library(scoringutils)
source("src/plot_functions.R")

model_info <-
  read_csv("model_table.csv") %>%
  distinct(model_design, .keep_all = TRUE) %>%
  select(-c(model_id, seed)) %>%
  left_join(
    tibble(file_path = dir_ls("results/posterior_predictive/")) %>%
      mutate(model_design = file_path %>%
        str_extract("(?<=model_design=)\\d+") %>%
        as.integer()) %>%
      select(-file_path)
  ) %>%
  distinct(
    constant_alpha,
    constant_IFR,
    constant_R0,
    double_IFR_0,
    half_alpha_0,
    half_R0_0,
    half_S_0,
    use_seroprev,
    use_tests
  ) %>%
  mutate(model_group = seq_len(n()))

pred_score_tbl <- map_dfr(dir_ls("results/prediction_score"), read_csv)

forecast_crps_plot <-
  pred_score_tbl %>%
  add_count(model) %>%
  filter(n > 8) %>%
  select(-n) %>%
  filter(
    target_type == "deaths",
    weeks_ahead %in% c(1, 4)
  ) %>%
  left_join(model_info, by = c("model" = "model_group")) %>%
  ggplot(aes(date, crps, color = use_tests, linetype = use_seroprev)) +
  facet_wrap(. ~ weeks_ahead,
    ncol = 1, scales = "free_y",
    labeller = as_labeller(
      ~ str_c(., " Week", if_else(. == 1, "", "s"), " Ahead Forecasts")
    )
  ) +
  geom_line() +
  geom_point() +
  theme_minimal_grid() +
  scale_x_date("Prediction Target Date") +
  scale_y_continuous("Continuous Ranked Probability Score (Lower is better)") +
  scale_color_discrete("Test Data",
    labels = ~ str_c(if_else(., "", "Not "), "Conditioned on Tests")
  ) +
  scale_linetype_discrete("Seroprev Data",
    labels = ~ str_c(if_else(., "", "Not "), "Using Seroprevalence Data")
  ) +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  guides(
    color = guide_legend(reverse = TRUE),
    linetype = guide_legend(reverse = TRUE)
  )


save_plot(
  filename = "/Users/damon/Documents/semi_parametric_COVID_19_OC_manuscript/figures/forecast_crps_plot.pdf",
  plot = forecast_crps_plot,
  ncol = 1,
  nrow = 2,
  base_asp = 3
)
