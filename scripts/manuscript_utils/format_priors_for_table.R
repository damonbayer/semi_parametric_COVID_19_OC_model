library(tidyverse)
library(scales)
expit <- function(x)  1/(1+exp(-x))
logit <- function(x) -log(1/x - 1)

format_for_table <- function(dist_type = "log", dist_mean, dist_sd, interval_width = 0.95) {
  tail_width <- (1 - interval_width) / 2 
  raw_quantiles <- qnorm(p = c(0.5, tail_width, 1 - tail_width), mean = dist_mean, sd = dist_sd)
  if (dist_type == "log") {
    trans_quantiles <- exp(raw_quantiles)
    distribuion_name <- "Log-Normal"
  } else if (dist_type == "logit") {
    trans_quantiles <- expit(raw_quantiles)
    distribuion_name <- "Logit-Normal"
  }
  rounded_dist_mean <- dist_mean %>% signif(3) %>% format(scientific = FALSE)
  rounded_dist_var <- dist_sd^2 %>% signif(3) %>% format(scientific = FALSE)
  rounded_quantiles <- trans_quantiles %>% signif(3) %>% format(scientific = FALSE)
  str_c(
    distribuion_name,
    "(",
    rounded_dist_mean,
    ", ",
    rounded_dist_var,
    ")",
    " & ",
    "\\makecell{", 
    rounded_quantiles[1], 
    " \\\\ (", 
    rounded_quantiles[2],
    ", ",
    rounded_quantiles[3],
    ")}") %>% 
    cat()
}

format_for_table(dist_type = "log", dist_mean =  1.35 - log(2), dist_sd = 0.11)
