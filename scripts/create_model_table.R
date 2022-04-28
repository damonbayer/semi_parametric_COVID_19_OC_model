library(tidyverse)
model_table <-
  bind_rows(
    # forecasts with and without seroprev
    crossing(
      max_t = 20:42,
      use_seroprev = c(T, F)
    ),
    # forecasts with and without tests
    crossing(
      max_t = 20:42,
      use_tests = c(T, F)
    ),
    # sensitivity with and without seroprevalence
    crossing(
      use_seroprev = c(T, F),
      half_S_0 = c(T, F)
    ),
    crossing(
      use_seroprev = c(T, F),
      half_R0_0 = c(T, F)
    ),
    crossing(
      use_seroprev = c(T, F),
      double_IFR_0 = c(T, F)
    ),
    crossing(
      use_seroprev = c(T, F),
      half_alpha_0 = c(T, F)
    ),
    # constant time-varying parameters
    tibble(constant_R0 = T),
    tibble(constant_alpha = T),
    tibble(constant_IFR = T)
  ) %>%
  # replace missing values with defualts
  replace_na(list(
    max_t = 42,
    use_seroprev = T,
    use_tests = T,
    half_S_0 = F,
    half_R0_0 = F,
    double_IFR_0 = F,
    half_alpha_0 = F,
    constant_R0 = F,
    constant_alpha = F,
    constant_IFR = F
  )) %>%
  distinct() %>%
  mutate(model_design = row_number(), .before = 1) %>%
  crossing(seed = 1:4) %>%
  mutate(model_id = row_number(), .before = 1) %>%
  relocate(sort(tidyselect::peek_vars())) %>%
  select(model_id, model_design, seed, everything())

write_csv(model_table, "model_table.csv")
