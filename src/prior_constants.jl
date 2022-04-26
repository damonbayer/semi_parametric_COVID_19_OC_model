## Cases
if use_tests
  σ_α_non_centered_mean = -3
  σ_α_non_centered_sd = 0.2
  α_init_non_centered_sd = 0.15
  if half_alpha_0
    α_init_non_centered_mean = 1.35 - log(2)
  else
    α_init_non_centered_mean = 1.35
  end
  else
    σ_ρ_cases_non_centered_mean = -2.5
    σ_ρ_cases_non_centered_sd = 0.2
    ρ_cases_init_non_centered_mean = -2.5
    ρ_cases_init_non_centered_sd = 0.2   
  end

## R0
if half_R0_0
  R₀_init_non_centered_mean = 0.0 - log(2)
else
  R₀_init_non_centered_mean = 0.0
end
R₀_init_non_centered_sd = 0.5
σ_R0_non_centered_mean = -3
σ_R0_non_centered_sd = 0.2

## IFR
if double_IFR_0
  IFR_init_non_centered_mean = -5.3 + log(2)
else
  IFR_init_non_centered_mean = -5.3
end
IFR_init_non_centered_sd = 0.4
σ_IFR_non_centered_mean = -2.5
σ_IFR_non_centered_sd = 0.13

## Intitial Conditions
if half_S_0
  S_SEI_non_centered_mean = 6 - log(2)
else
  S_SEI_non_centered_mean = 6
end
S_SEI_non_centered_sd = 0.5
I_EI_non_centered_mean = 0.6
I_EI_non_centered_sd = 0.03

## Dwell Times
# dur_latent_non_centered_mean = log(3 / 7)
# dur_latent_non_centered_sd = 0.25
# https://academic.oup.com/cid/advance-article/doi/10.1093/cid/ciab746/6359063
dur_latent_non_centered_mean = -0.25
dur_latent_non_centered_sd = 0.1

# dur_infectious_non_centered_mean = log(5 / 7)
# dur_infectious_non_centered_sd = 0.2
# https://bmjopen.bmj.com/content/10/8/e039856.abstract
dur_infectious_non_centered_mean = 0.15
dur_infectious_non_centered_sd = 0.1
ρ_death_non_centered_mean = 2.3
ρ_death_non_centered_sd = 0.2

## Overdispersion
include(projectdir("src/prior_constants_overdispersion.jl"))
