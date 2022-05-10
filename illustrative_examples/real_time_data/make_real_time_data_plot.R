library(tidyverse)
library(lubridate)
library(here)
library(data.table)
library(dtplyr)
library(scales)
library(cowplot)

negative_line_list_path_new <- "/Users/damon/Documents/uci_covid_modeling2/data/from_OCHCA/All PCR tests updated 05.10.22.csv"
negative_line_list_path_old <- "/Users/damon/Downloads/All PCR tests updated 1.6.22.csv"

process_neg_line_list <- function(negative_line_list_path) {
  neg_line_list <-
    fread(negative_line_list_path,
          # select = c("unique_num", "Specimen.Collected.Date", "TestResult", "Zip"),
          select = c("Specimen.Collected.Date", "TestResult", "Zip"),
          na.strings = "") %>%
    # select(id = unique_num, date = Specimen.Collected.Date, zip = Zip, test_result = TestResult) %>%
    select(date = Specimen.Collected.Date, zip = Zip, test_result = TestResult) %>%
    # filter(!(is.na(id) | is.na(test_result))) %>%
    filter(!is.na(test_result)) %>%
    mutate(test_result = fct_collapse(str_to_lower(test_result),
                                      negative = "negative",
                                      positive = "positive",
                                      other = c("invalid", "inconclusive"))) %>%
    filter(date >= lubridate::ymd("2020-01-01")) %>%
    # group_by(id) %>%
    arrange(date) %>%
    ungroup()

  neg_line_list_filtered <- neg_line_list

  oc_data <-
    neg_line_list_filtered %>%
    count(date, test_result) %>%
    pivot_wider(names_from = test_result, values_from = n, values_fill = 0) %>%
    as_tibble() %>%
    replace(is.na(.), 0) %>%
    mutate(cases = positive, tests = negative + positive + other) %>%
    select(date, cases, tests) %>%
    mutate(date = as_date(date))
  oc_data
}

lump_data <- function(oc_data,
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
                     .groups = "drop") %>%
    mutate(time = as.numeric((end_date - min(start_date) + 1) / 7)) %>%
    dplyr::select(time, everything(), -lump) %>%
    arrange(start_date)

  if((head(lumped_oc_data$start_date, 1) != first_day) | (tail(lumped_oc_data$end_date, 1) != last_day)) {
    warning("Lumped data has different range than input")
  }

  lumped_oc_data
}

new_dat <- process_neg_line_list(negative_line_list_path_new)
old_dat <- process_neg_line_list(negative_line_list_path_old)

last_day <- "2022-01-05"
first_day <- as.character(ymd(last_day) - 28)

tmp <-
  bind_rows(
  lump_data(new_dat, first_day, last_day, time_interval_in_days = 7) %>%
    mutate(reported_on = mdy("05.10.22"),
           lump = 7),
  lump_data(old_dat, first_day, last_day, time_interval_in_days = 7) %>%
    mutate(reported_on = mdy("1.06.21"),
           lump = 7),
  lump_data(new_dat, first_day, last_day, time_interval_in_days = 3) %>%
    mutate(reported_on = mdy("05.10.22"),
           lump = 3),
  lump_data(old_dat, first_day, last_day, time_interval_in_days = 3) %>%
    mutate(reported_on = mdy("1.06.21."),
           lump = 3)) %>%
  mutate(test_positivity = cases / tests) %>%
  select(date = end_date, cases, tests, test_positivity, reported_on, binning = lump) %>%
  pivot_longer(c(cases, tests, test_positivity)) %>%
  mutate(reported_on = format(reported_on, "%b %d, %y"))

tmp_plot <-
  tmp %>%
  mutate(binning = str_c(binning, " days"),
         data = name %>%
           fct_relevel(c("tests", "cases", "test_positivity")) %>%
           fct_relabel(~str_replace(., "_",  " ") %>%
                         str_to_title())) %>%
  rename_with(str_to_title, c(binning, data)) %>%
  ggplot(aes(date, value, color = reported_on)) +
  facet_grid(Data ~ Binning, scales = "free_y", labeller = label_both) +
  geom_line() +
  geom_point() +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "bottom") +
  scale_x_date(name = "Date", date_labels = "%b %d '%y") +
  scale_y_continuous(NULL, labels = comma) +
  scale_color_discrete("Date Reported")



save_plot_target_asp <- function (filename, plot, ncol = 1, nrow = 1, base_height = 3.71,
                                  base_asp = 1.618, base_width = NULL) {
  cowplot::save_plot(filename, plot, ncol = ncol, nrow = nrow, base_height = base_height,
                      base_asp = base_asp * nrow / ncol, base_width = base_width)
}


save_plot_target_asp(filename = "figures/advancement_slides/real_time_data.pdf",
                     plot = tmp_plot,
                     ncol = 2,
                     nrow = 3,
                     base_asp = 16/9,
                     base_height = 2)
