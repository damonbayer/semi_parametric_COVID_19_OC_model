library(gridExtra)
library(tidyverse)
library(tidybayes)
library(fs)
library(scales)
library(latex2exp)
library(cowplot)

dat <- 
  read_csv("data/oc_data.csv") %>% 
  mutate(index = 1:n()) %>% 
  rename(date = end_date) %>% 
  select(-start_date)

dat_tidy <-
  dat %>% 
  select(date, cases, tests, deaths) %>% 
  mutate(test_positivity = cases / tests) %>% 
  select(-tests) %>% 
  pivot_longer(-date)

model_table <-
  read_csv("model_table.csv") %>% 
  select(-model_id, -seed) %>% 
  distinct()

all_pp <- 
  tibble(full_path = dir_ls("results/tidy_posterior_predictive")) %>% 
  mutate(model_design = full_path %>% 
           str_extract("(?<=model_design=)\\d+") %>% 
           as.integer()) %>% 
  left_join(model_table) %>% 
  mutate(pp_data = map(full_path, read_csv)) %>% 
  select(-full_path) %>% 
  arrange(model_design)


make_pp_plot_comparison <- function(plot_max_t, n_forecast_times = 4) {
  pp_tmp <- 
    all_pp %>% 
    filter(max_t == plot_max_t,
           constant_alpha == F,
           constant_IFR == F,
           constant_R0 == F,
           double_IFR_0 == F,
           half_alpha_0 == F,
           half_R0_0 == F,
           half_S_0 == F) %>% 
    mutate(model_description = str_c("use_seroprev=", use_seroprev, "\n", "use_tests=", use_tests)) %>% 
    select(model_description, max_t, pp_data) %>% 
    unnest(pp_data) %>% 
    filter(name %in% c("deaths", "test_positivity"))
  
  first_forecast_date <-
    pp_tmp$date %>%
    unique() %>% 
    sort() %>% 
    tail(n_forecast_times) %>% 
    head(1)
  
  dat_tidy_tmp <- 
    dat_tidy %>% 
    filter(name %in% c("deaths", "test_positivity")) %>% 
    mutate(forecast = date >= first_forecast_date)
  
  ggplot(mapping = aes(date, value)) +
    facet_grid(name ~ model_description, scales = "free_y") +
    geom_lineribbon(data = pp_tmp,
                    mapping = aes(ymin = .lower, ymax = .upper)) +
    geom_point(data = dat_tidy_tmp,
               mapping = aes(shape = forecast, color = forecast)) +
    scale_fill_brewer(name = "Credible Interval Width",
                      labels = ~percent(as.numeric(.))) +
    theme_minimal_grid() +
    theme(legend.position = "bottom",
          plot.title = element_text(size=12)) +
    guides(fill = guide_legend(reverse = TRUE)) +
    ggtitle(str_c("max_t=", plot_max_t))
} 


all_pp_comparison_plots <- map(sort(unique(all_pp$max_t)), make_pp_plot_comparison)

ggsave2(filename = "figures/all_pp_comparison_plots.pdf",
        plot = marrangeGrob(all_pp_comparison_plots, nrow=1, ncol=1),
        width = 8,
        height = 8)
