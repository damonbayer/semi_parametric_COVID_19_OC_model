gq <- read_csv("multiple_map_gq_chains.csv")


lp <- c(-587.4482515100262,
-587.680340164783,
-587.7895387027363,
-588.2706537244026,
-588.6108127679697,
-589.0281953687298,
-589.7019716383043,
-589.8378180392906,
-590.213309066774,
-590.5094751524744)



tmp <- 
  gq %>% 
  pivot_longer(-c(iteration, chain), names_to = "name_raw") %>% 
  separate(name_raw, c("name", "time"), sep = "\\[", remove = T) %>% 
  drop_na() %>% 
  mutate(time = as.integer(str_sub(time, end = -2))) %>% 
  mutate(name = str_remove(name, "data_new_")) %>% 
  mutate(lp = lp[iteration],
         seed = as_factor(iteration))

tmp %>% 
ggplot(aes(time, value, color = lp, group = lp)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_color_viridis_c(option = "D", direction = -1) +
  ggtitle("Time Varying Generated Quantities for Top 10 MAP Estimates") +
  theme(legend.position = c(0.5, 0.1))

                           