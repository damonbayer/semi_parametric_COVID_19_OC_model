library(tidyverse)

model_table <- 
  bind_rows(
    tibble(max_t = 20:42),
    tibble(use_seroprev = F),
    tibble(use_tests = F),
    tibble(half_S_0 = T),
    tibble(half_R0_0 = T),
    tibble(double_IFR_0 = T),
    tibble(half_alpha_0 = T),
    tibble(constant_R0 = T),
    tibble(constant_alpha = T),
    tibble(constant_IFR = T)
    ) %>%
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
    constant_IFR = F)) %>%
  mutate(model_design = row_number()) %>% 
  crossing(seed = 1:4) %>%
  mutate(model_id = row_number()) %>%
  select(model_id, model_design, seed, everything())

write_csv(model_table, "model_table.csv")

model_table %>% 
  distinct(model_design, .keep_all = T) %>% 
  View()

crossing(max_t = 20:42,
         use_seroprev = c(T, F),
         use_tests = c(T, F),
         # half_S_0 = c(T, F),
         # half_R0_0 = c(T, F),
         # double_IFR_0 = c(T, F),
         # half_alpha_0 = c(T, F),
         constant_R0 = c(T, F),
         constant_alpha = c(T, F),
         constant_IFR = c(T, F))

