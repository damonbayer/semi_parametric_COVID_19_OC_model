library(gridExtra)
library(tidyverse)
library(tidybayes)
library(fs)
library(scales)
library(latex2exp)
library(cowplot)

model_table <- read_csv("model_table.csv")
n_forecast_times <- 4

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


make_pp_plot <- function(path, n_forecast_times = 4) {
  file_description <- 
    path %>% 
    path_file() %>% 
    str_remove("tidy_posterior_predictive") %>% 
    str_extract_all("(?<=_)\\w+=[:alnum:]+(?=_|\\.)") %>% 
    unlist() %>%
    str_c(collapse = "\n")
  
  pp <- read_csv(path)  
  
  first_forecast_date <-
    pp$date %>%
    unique() %>% 
    sort() %>% 
    tail(n_forecast_times) %>% 
    head(1)
  
  ggplot(mapping = aes(date, value)) +
    facet_wrap(. ~ name, scales = "free_y") +
    geom_lineribbon(data = pp %>% 
                      filter(name %in% c("deaths", "test_positivity")),
                    mapping = aes(ymin = .lower, ymax = .upper)) +
    geom_point(data = dat_tidy %>% 
                 filter(name %in% c("deaths", "test_positivity")) %>% 
                 mutate(forecast = date >= first_forecast_date),
               mapping = aes(shape = forecast, color = forecast)) +
    scale_fill_brewer(name = "Credible Interval Width",
                      labels = ~percent(as.numeric(.))) +
    theme_minimal_grid() +
    theme(legend.position = "bottom",
          plot.title = element_text(size=12)) +
    guides(fill = guide_legend(reverse = TRUE)) +
    ggtitle(file_description)
}


all_pp_plots <- map(dir_ls("results/tidy_posterior_predictive"), make_pp_plot)

ggsave2(filename = "figures/all_pp_plots.pdf",
        plot = marrangeGrob(all_pp_plots, nrow=1, ncol=1),
        width = 8,
        height = 8)
