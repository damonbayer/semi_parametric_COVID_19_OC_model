library(tidyverse)
library(lubridate)
oc_data <- read_csv("~/Documents/uci_covid_modeling2/data/oc_data.csv")
death_delay_ecdf <- read_rds("~/Documents/uci_covid_modeling2/data/death_delay_ecdf.rds")
first_day <- ymd("2020-03-30")
last_day <- ymd("2021-02-28")

lump_oc_data <- function(oc_data,
                         time_interval_in_days,
                         first_day = "0000-01-01",
                         last_day = "9999-12-31") {
  
  lumped_oc_data <- oc_data %>%
    filter(date >= lubridate::ymd(first_day),
           date <= lubridate::ymd(last_day)) %>%
    group_by(lump = as.integer(floor((max(date) - date) / time_interval_in_days))) %>%
    filter(n() == time_interval_in_days) %>%
    dplyr::summarize(start_date = min(date),
                     end_date = max(date),
                     cases = sum(cases),
                     tests = sum(tests),
                     deaths = sum(deaths),
                     .groups = "drop") %>%
    mutate(time = as.numeric((end_date - min(start_date) + 1) / 7)) %>%
    dplyr::select(time, everything(), -lump) %>%
    arrange(start_date)
  
  if((head(lumped_oc_data$start_date, 1) != first_day) | (tail(lumped_oc_data$end_date, 1) != last_day)) {
    warning("Lumped data has different range than input")
  }
  
  lumped_oc_data
}

dat <- oc_data %>%
  lump_oc_data(time_interval_in_days = 7,
               first_day,
               last_day) %>%
  mutate(prop_deaths_reported = death_delay_ecdf(as.numeric(max(oc_data$date) - end_date)))

dat_seroprev <- tibble(start_date = c(ymd("2020-07-10"), ymd("2021-01-30")),
                       end_date = c(ymd("2020-08-16"),  ymd("2021-02-20")),
                       seroprev_cases = c(343, 224),
                       seroprev_tests = c(2979, 495)) %>%
  mutate(time = round(as.numeric((end_date - min(dat$start_date) + 1) / 7))) %>%
  select(time, everything()) %>%
  slice(1)


write_csv(dat, "data/oc_data.csv")
write_csv(dat_seroprev, "data/oc_seroprev_data.csv")
