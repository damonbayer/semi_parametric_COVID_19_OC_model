library(tidyverse)
source("src/plot_functions.R")

all_metrics <-
  read_csv("results/simulation/oc_like/simulated_rt_comparison_all_metrics.csv") %>%
  mutate(method = method |> 
           fct_inorder() |> 
           fct_recode(`True Model` = "Full Model")) %>% 
  filter(method != "Rt-estim-gamma")

envelope_plot <-
  all_metrics %>%
  ggplot(aes(method, mean_env)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  labs(title = "Envelope", x = "Method", y = "Envelope") +
  theme_minimal_grid()

mciw_plot <-
  all_metrics %>%
  ggplot(aes(method, MCIW)) +
  geom_boxplot() +
  labs(title = "MCIW", x = "Method", y = "MCIW") +
  theme_minimal_grid()

abs_dev_plot <-
  all_metrics %>%
  ggplot(aes(method, mean_dev)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Absolute Deviation", x = "Method", y = "Deviation") +
  theme_minimal_grid()

masv_plot <-
  all_metrics %>%
  ggplot(aes(method, MASV)) +
  geom_boxplot() +
  geom_hline(yintercept = unique(all_metrics$true_MASV), linetype = "dashed") +
  labs(title = "MASV", x = "Method", y = "MASV") +
  theme_minimal_grid()

rt_comparison_metrics_plot <- plot_grid(envelope_plot, mciw_plot, abs_dev_plot, masv_plot, align = "hv", nrow = 2, ncol = 2)

save_plot(filename = path(figures_dir, "rt_comparison_metrics_plot", ext = "pdf"),
          plot = rt_comparison_metrics_plot, ncol = 2, nrow = 2)
