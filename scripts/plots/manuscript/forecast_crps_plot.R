library(tidyverse)
library(fs)
library(scoringutils)
source("src/plot_functions.R")

model_info <-
  read_csv("model_table.csv") %>%
  select(-model_id, -seed) %>%
  distinct()

models_to_keep <-
  model_info %>%
  count(
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
  filter(n > 1) %>%
  mutate(model_group = 1:n()) %>%
  left_join(model_info) %>%
  rename(model = model_design)

pred_score_tbl <- map_dfr(dir_ls("results/posterior_predictive_score"), read_csv) %>%
  arrange(model)

forecast_crps_plot <-
  pred_score_tbl %>%
  filter(model %in% models_to_keep$model) %>%
  filter(
    target_type == "deaths",
    weeks_ahead %in% c(1, 4)
  ) %>%
  left_join(models_to_keep) %>%
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
  ) +
  ggtitle("Continuous Ranked Probability Score Comparison")


save_plot(
  filename = "~/Documents/semi_parametric_COVID_19_OC_manuscript/figures/forecast_crps_plot.pdf",
  plot = forecast_crps_plot,
  ncol = 1,
  nrow = 2,
  base_asp = 2
)
