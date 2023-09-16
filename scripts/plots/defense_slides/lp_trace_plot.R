# need an lp plot
library(tidyverse)
source("src/plot_functions.R")

target_model_design <- 46

posterior_lp <- 
  dir_ls("results/posterior_lp") %>% 
  enframe(name = NULL) %>% 
  filter(str_detect(value, str_c("model_design=", target_model_design, "_"))) %>% 
  pull(value) %>% 
  read_csv()

lp_trace_plot <- 
  posterior_lp %>% 
  mutate(.chain = factor(.chain)) %>% 
  ggplot(aes(.iteration, lp, color = .chain)) +
  geom_line() +
  labs(title = "Log-Posterior Traceplot",
       x = "Iteration",
       y = "Log-Posterior",
       color = "Chain")

save_plot_target_asp(filename = path(defense_figures_dir, "lp_trace_plot", ext = "pdf"),
          plot = lp_trace_plot)

