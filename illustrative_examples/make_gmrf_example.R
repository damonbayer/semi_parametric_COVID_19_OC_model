# Illustrate relationship between GMRF prior and corresponding time-varying function
library(tidyverse)
library(latex2exp)
source("src/plot_functions.R")
set.seed(4)

expit <- function(x)  1/(1+exp(-x))
logit <- function(x) -log(1/x - 1)

ifr_bar_tbl <-
  tibble(ifr_bar = cumsum(
    c(logit(0.05),
      rnorm(5, sd = 0.3)))) %>%
  mutate(l = row_number())

ifr_t_tbl <- ifr_bar_tbl %>%
  transmute(ifr_t = expit(ifr_bar),
            t = (l - 1) * 7)

gmrf_example <-
  cowplot::plot_grid(
    ggplot(ifr_bar_tbl, aes(l, ifr_bar)) +
      geom_point() +
      scale_y_continuous(name = TeX("$\\tilde{\\eta}_l$")) +
      scale_x_continuous(name = "l (index)", breaks = 1:7) +
      cowplot::theme_minimal_grid(),
    ggplot(ifr_t_tbl, aes(t, ifr_t)) +
      geom_step() +
      geom_point() +
      scale_y_continuous(name = TeX("$\\eta(t)$"), labels = scales::percent) +
      scale_x_continuous(name = "t (days)", breaks = seq(0, 35, 7)) +
      cowplot::theme_minimal_grid(),
    ncol = 2, align = "hv")

save_plot_target_asp(filename = "figures/advancement_slides/gmrf_example.pdf",
          plot = gmrf_example,
          ncol = 2, nrow = 1,
          base_asp = 32/9)
