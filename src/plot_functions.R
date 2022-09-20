library(ggplot2)
library(cowplot)
library(scales)
library(latex2exp)
library(fs)
theme_set(theme_minimal_grid())
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
    "dur_infectious_days" = "$1/\\nu$",
    "dur_latent_days" = "$1/\\gamma$",
    "E" = "E",
    "I" = "I",
    "I_EI" = "Initial $\\frac{I}{E + I}$",
    "IFR_t" = "$\\eta$",
    "IFR" = "$\\eta$",
    "R" = "R",
    "R₀_t" = "$R_0$",
    "R₀" = "$R_0$",
    "Rₜ_t" = "$R_t$",
    "S" = "S",
    "S_SEI" = "Initial $\\frac{S}{S + E + I}$",
    "seroprev_mean" = "Seroprev Mean",
    "α_t" = "$\\alpha$",
    "α" = "$\\alpha$",
    "β_t" = "$\\beta$",
    "β" = "$\\beta$",
    "ρ_cases_t" = "$\\rho^Y$",
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
