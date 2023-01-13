library(ggplot2)
library(cowplot)
library(scales)
library(latex2exp)
library(fs)
library(ggdist)
library(ggblend)
theme_set(theme_minimal_grid())
popsize <- 3175692
save_plot_target_asp <- function (filename, plot, ncol = 1, nrow = 1, base_height = 3.71,
                                  base_asp = 1730/650, base_width = NULL) {
  cowplot::save_plot(filename, plot, ncol = ncol, nrow = nrow, base_height = base_height,
                     base_asp = base_asp * nrow / ncol, base_width = base_width)
}

figures_dir <- path("~/Documents/semi_parametric_COVID_19_OC_manuscript/figures")

my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme(legend.position = "bottom"))

brewer_line_color <- "#08519c"

my_labeller <- 
  c("cases_bb_mean" = "Cases BB Mean",
    "cases_mean" = "Cases Mean",
    "cases_nb_mean" = "Cases NB Mean",
    "D" = "D",
    "deaths_mean" = "Deaths Mean",
    # "dur_infectious_days" = "$Infectious\\;Period\\;(Days)$",
    # "dur_latent_days" = "$Latent\\;Period\\;(Days)$",
    "dur_infectious_days" = "$7/\\nu$",
    "dur_latent_days" = "$7/\\gamma$",
    "E" = "E",
    "I" = "I",
    # "I_EI" = "Initial $\\frac{I}{E + I}$",
    "I_EI" = "$\\tilde{I}_{0}$",
    "IFR_t" = "$\\eta(t)$",
    "IFR" = "$\\eta$",
    "R" = "R",
    "R₀_t" = "$R_0(t)$",
    "R₀" = "$R_0$",
    "Rₜ_t" = "$R_t(t)$",
    "S" = "S",
    # "S_SEI" = "Initial $\\frac{S}{S + E + I}$",
    "S_SEI" = "$S_0$",
    "seroprev_mean" = "Seroprev Mean",
    "α_t" = "$\\alpha(t)$",
    "α" = "$\\alpha$",
    "β_t" = "$\\beta(t)$",
    "β" = "$\\beta$",
    "ρ_cases_t" = "$\\rho^Y(t)$",
    "ρ_death" = "$\\rho^D$",
    "σ_IFR" = "$\\sigma_\\eta$",
    "σ_R0" = "$\\sigma_{R_0}$",
    "σ_α" = "$\\sigma_\\alpha$",
    "σ_ρ_cases" = "$\\sigma^2_{\\rho^Y}",
    "ϕ_cases_bb" = "$\\phi_C$",
    "ϕ_cases_nb" = "$\\phi_Y$",
    "ϕ_deaths" = "$\\phi_D$") %>%
  TeX(output = "expression")

my_labeller_fn <- as_labeller(function(string) my_labeller[string], label_parsed)
