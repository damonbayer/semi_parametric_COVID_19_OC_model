coverage_summary <- read_csv("results/coverage_summary.csv")

coverage_summary %>% 
  filter(str_detect(name, "\\[\\d+\\]", negate = T)) %>% 
  drop_na() %>% 
  filter(.width == 0.8) %>% 
  arrange(desc(covered))

true_parameters$name

true_parameters %>% 
  filter(str_detect(name, "σ"))

coverage_summary %>% 
  mutate(time = str_extract(name, "(?<=\\[)\\d+(?=\\])") %>% as.numeric()) %>% 
  mutate(name = str_extract(name, "^.+(?=\\[)")) %>% 
  drop_na() %>% 
  ggplot(aes(time, covered)) +
  # facet_wrap(.width ~ name) +
  facet_wrap(name ~ .width, labeller = label_both) +
  geom_line()

coverage_summary %>% 
  drop_na() %>% 
  group_by(.width) %>%
  summarize(mean(covered < .width))


coverage_summary %>% 
  ggplot(aes(covered)) +
  facet_wrap(. ~  .width) +
  geom_histogram(binwidth = 0.05)
  
coverage_summary %>% 
  filter(.width == 0.8) %>% 
  arrange(covered) %>% 
  ggplot(aes(covered)) +
  geom_histogram()


tmp <- 
  tibble(full_path = dir_ls("results/simulated_posterior_samples_csv")) %>% 
  mutate(tidy_summary = map(full_path,
                            ~read_csv(.) %>% 
                              select(-iteration, -chain) %>% 
                              pivot_longer(everything()) %>% 
                              group_by(name) %>% 
                              median_qi(.width = c(0.5, 0.8, 0.95)))) %>% 
  unnest(tidy_summary) %>% 
  select(-full_path) %>% 
  left_join(true_parameters) %>% 
  mutate(lower_error = .lower > true_value,
         upper_error = true_value > .upper) %>% 
  mutate(covered = !lower_error & !upper_error)

tmp %>% 
  filter(name == "ϕ_cases_non_centered")
