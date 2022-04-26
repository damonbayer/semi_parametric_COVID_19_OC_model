library(tidyverse)
library(stemr)
library(tidybayes)

n_sims <- 2000
n_t <- 40

init_mean <- -2.5
init_sd <- 0.2
step_sd <- 0.1

generate_cdr_t <- function(x) expit(rnorm(1, init_mean, init_sd) + c(0, cumsum(rnorm(n_t - 1, 0, step_sd))))

rerun(n_sims, generate_cdr_t()) %>% 
  enframe(name = "iter") %>% 
  unnest(value) %>% 
  group_by(iter) %>% 
  mutate(t = row_number()) %>% 
  group_by(t) %>% 
  select(-iter) %>% 
  mutate(value = 1/value) %>% 
  median_qi(.width = c(0.5, 0.8, 0.95)) %>% 
  ggplot(aes(t, value, ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = c(3, 5, 10, 20))

