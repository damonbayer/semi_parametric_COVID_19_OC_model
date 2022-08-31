library(tidyverse)
library(fs)
library(colorblindr)
library(cowplot)
library(tidybayes)

model_table <- read_csv("model_table.csv") %>% 
  distinct(model_design, .keep_all = T) %>% 
  left_join(read_csv("data/oc_data.csv", col_select = c(max_t = time, date = end_date))) %>% 
  select(-seed, -model_id)


mean_log_lik_tbl <- 
  tibble(file_name = dir_ls("results/prediction_score")) %>% 
  mutate(model_design = file_name %>%
           str_extract("(?<=model_design=)\\d+") %>% 
           as.integer()) %>% 
  left_join(model_table) %>% 
  mutate(prediction_score = map2(file_name, max_t,
                                 ~read_csv(.x) %>% 
                                   drop_na(date) %>% 
                                   mutate(time = 1:n()) %>% 
                                   mutate(weeks_ahead = as.integer(time - .y)) %>% 
                                   filter(weeks_ahead > 0) %>% 
                                   select(-time))) %>% 
  select(-date) %>%
  arrange(model_design) %>% 
  unnest(prediction_score) %>% 
  pivot_longer(cols = starts_with("ll_"),
               names_prefix = "ll_",
               values_to = "log_lik") %>% 
  drop_na() %>% 
  group_by(model_design, weeks_ahead, name) %>% 
  summarize(mean_log_lik = mean(log_lik)) %>% 
  left_join(model_table)


mean_log_lik_tbl %>% 
  filter(max_t == max(max_t)) %>% 
  ggplot(aes(x = mean_log_lik, fill = use_tests)) +
  facet_grid(weeks_ahead ~ name, scales = "free_x") +
  geom_dotsinterval(alpha = 0.5) +
  cowplot::theme_cowplot()

mean_log_lik_tbl %>% 
  filter(max_t == max(max_t)) %>% 
  ggplot(aes(x = mean_log_lik, fill = use_seroprev)) +
  facet_grid(weeks_ahead ~ name, scales = "free_x") +
  geom_density(alpha = 0.5) +
  cowplot::theme_cowplot()


tmp %>% 
  ggplot(aes(x = mean_log_lik, y = ".", label = model_design)) +
  facet_grid(weeks_ahead ~ name, scales = "free_x") +
  geom_text() +
  cowplot::theme_cowplot()
  
tmp
tmp %>% select(-file_name) %>% unnest(prediction_score)


tmp %>% 
  filter(constant_alpha == F,
         constant_IFR == F,
         constant_R0 == F,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F) %>% 
  ggplot(aes(x = max_t, y = mean_log_lik, color = use_tests, shape = use_seroprev)) +
  facet_grid(weeks_ahead ~ name, scale = "free_y") +
  geom_point() +
  geom_line()



# -------------------------------------------------------------------------
pred_score_plot <- 
  mean_log_lik_tbl %>% 
  left_join(model_table) %>% 
  filter(constant_alpha == F,
         constant_IFR == F,
         constant_R0 == F,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F) %>% 
  mutate(category_name = str_c(if_else(use_tests, "Condition on Tests", "Not Condition on Tests"),
                               "\n",
                               if_else(use_seroprev, "Use Seroprevalence Data", "Not Use Seroprevalence Data"),
                               sep = " ")) %>% 
  ggplot(aes(date, mean_log_lik, color = category_name)) +
  facet_grid(weeks_ahead ~ name, scales = "free",
             labeller = labeller(name = ~str_to_title(.),
                                 weeks_ahead = ~str_c(., "Weeks Ahead", sep = " "))) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous("Mean Log-Likelihood") +
  scale_x_date("Last Fit Date") +
  scale_color_OkabeIto(name = "Model Choices") +
  theme(legend.position = "bottom") +
  ggtitle("Prediction Score Comparison")

save_plot(filename = "~/Desktop/pred_score_plot.pdf", plot = pred_score_plot, ncol = 2, nrow = 4)


pred_score_no_sero_plot <-
  mean_log_lik_tbl %>% 
  left_join(model_table) %>% 
  filter(constant_alpha == F,
         constant_IFR == F,
         constant_R0 == F,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F,
         use_seroprev == T) %>% 
  mutate(category_name = str_c(if_else(use_tests, "Conditioned on Tests", "Not Conditioned on Tests"))) %>% 
  ggplot(aes(date, mean_log_lik, color = category_name)) +
  facet_grid(weeks_ahead ~ name, scales = "free",
             labeller = labeller(name = ~str_to_title(.),
                                 weeks_ahead = ~str_c(., "Weeks Ahead", sep = " "))) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous("Mean Log-Likelihood") +
  scale_x_date("Last Fit Date") +
  scale_color_OkabeIto(name = "Model Choices") +
  theme(legend.position = "bottom") +
  ggtitle("Prediction Score Comparison", subtitle = "Not using Seroprevalence Data")

save_plot(filename = "~/Desktop/pred_score_no_sero_plot.pdf", plot = pred_score_no_sero_plot, ncol = 2, nrow = 4)


pred_score_sero_plot <-
  mean_log_lik_tbl %>% 
  left_join(model_table) %>% 
  filter(constant_alpha == F,
         constant_IFR == F,
         constant_R0 == F,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F,
         use_seroprev == T) %>% 
  mutate(category_name = str_c(if_else(use_tests, "Conditioned on Tests", "Not Conditioned on Tests"))) %>% 
  ggplot(aes(date, mean_log_lik, color = category_name)) +
  facet_grid(weeks_ahead ~ name, scales = "free",
             labeller = labeller(name = ~str_to_title(.),
                                 weeks_ahead = ~str_c(., "Weeks Ahead", sep = " "))) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous("Mean Log-Likelihood") +
  scale_x_date("Last Fit Date") +
  scale_color_OkabeIto(name = "Model Choices") +
  theme(legend.position = "bottom") +
  ggtitle("Prediction Score Comparison", subtitle = "Using Seroprevalence Data")

save_plot(filename = "~/Desktop/pred_score_sero_plot.pdf", plot = pred_score_sero_plot, ncol = 2, nrow = 4)

library(gt)
log_lik_tbl <- 
  tibble(file_name = dir_ls("results/prediction_score")) %>% 
  mutate(model_design = file_name %>%
           str_extract("(?<=model_design=)\\d+") %>% 
           as.integer()) %>% 
  left_join(model_table) %>% 
  mutate(prediction_score = map2(file_name, max_t,
                                 ~read_csv(.x) %>% 
                                   drop_na(date) %>% 
                                   mutate(time = 1:n()) %>% 
                                   mutate(weeks_ahead = as.integer(time - .y)) %>% 
                                   filter(weeks_ahead > 0) %>% 
                                   select(-time))) %>% 
  rename(max_t_date = date) %>%
  arrange(model_design) %>% 
  unnest(prediction_score) %>% 
  pivot_longer(cols = starts_with("ll_"),
               names_prefix = "ll_",
               values_to = "log_lik") %>% 
  drop_na()


mean_log_lik <- 
  log_lik_tbl %>% 
  select(-file_name) %>% 
  filter(constant_alpha == F,
         constant_IFR == F,
         constant_R0 == F,
         double_IFR_0 == F,
         half_alpha_0 == F,
         half_R0_0 == F,
         half_S_0 == F,
         name == "deaths") %>% 
  mutate(post_seroprev = date >= lubridate::ymd("2020-08-16")) %>% 
  group_by(post_seroprev, use_seroprev, use_tests, weeks_ahead) %>% 
  summarize(mean_log_lik = mean(log_lik), .groups = "drop") %>% 
  filter(weeks_ahead >= 2)

make_mean_log_lik_tbl <- function(post_seroprev_arg = T, caption = "change the caption") {
  mean_log_lik %>% 
    filter(post_seroprev == post_seroprev_arg) %>% 
    select(-post_seroprev) %>% 
    mutate(across(starts_with("use"), ~if_else(., "Yes", "No"))) %>% 
    mutate(weeks_ahead = str_c(weeks_ahead, "Weeks Ahead", sep = " ")) %>% 
    pivot_wider(names_from = weeks_ahead, values_from = mean_log_lik) %>% 
    gt(caption = caption) %>% 
    gt::tab_spanner("Mean Log-Likelihood of Deaths Forecast", columns = ends_with("Weeks Ahead")) %>%
    gt::tab_spanner("Model", columns = starts_with("use")) %>% 
    fmt_number(columns = ends_with("Weeks Ahead"), decimals = 2) %>% 
    cols_label(
      use_seroprev = "Use Seroprev Data",
      use_tests = "Use Tests")
}


gt_latex_dependencies() %>%
  as.character() %>%
  write_lines("~/Documents/semi_parametric_COVID_19_OC_manuscript/tables/gt_latex_dependencies.tex")

make_mean_log_lik_tbl(post_seroprev_arg = T,
                      caption = "Log-likelihood of deaths forecasts for pre-seroprevalence models") %>% 
  gtsave(filename = "~/Documents/semi_parametric_COVID_19_OC_manuscript/tables/log_lik_deaths_pre_sero_table.tex")

make_mean_log_lik_tbl(post_seroprev_arg = F,
                      caption = "Log-likelihood of deaths forecasts for post-seroprevalence models") %>% 
  gtsave(filename = "~/Documents/semi_parametric_COVID_19_OC_manuscript/tables/log_lik_deaths_post_sero_table.tex")

make_mean_log_lik_tbl(T) %>% 
  gtsave(filename = "my_table.tex")
make_mean_log_lik_tbl(F)
