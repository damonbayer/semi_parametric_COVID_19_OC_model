expit <- function(x)  1/(1+exp(-x))
logit <- function(x) -log(1/x - 1)
my_formatter <- scales::number_format(accuracy = 0.0001)


signif(6 - log(2), 3)
format_for_table <- function(x) {
  x_trans <- signif(x, 3)
  str_c("\\makecell{", 
        x_trans[1], 
        " \\\\ (", 
        x_trans[2],
        ", ",
        x_trans[3],
        ")}")
}




format_for_table(expit(qnorm(p = c(0.5, 0.025, 0.975), mean = 6 - log(2), sd = 0.5))) |> cat()
  
n <- 10000

tibble(alpha = rnorm(n = n, mean = 1.35, sd = 0.1),
         R0 = rnorm(n, mean = 0, sd = 0.25),
         IFR = rnorm(n, mean = -5.3, sd = 0.2)) |> 
  pivot_longer(everything(), names_to = "parameter", values_to = "y") |> 
  mutate(y = if_else(parameter == "IFR", expit(y), exp(y))) |> 
  ggplot(aes(y)) +
  facet_wrap(~parameter, scales = "free")+
  stat_halfeye(normalize = "panels")
