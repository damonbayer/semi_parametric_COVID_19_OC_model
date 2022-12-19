## Cases
if use_tests
  const σ_α_non_centered_mean = -2.7
  const σ_α_non_centered_sd = 0.15
  const α_init_non_centered_sd = 0.11
  if half_alpha_0
    const α_init_non_centered_mean = 1.35 - log(2)
  else
    const α_init_non_centered_mean = 1.35
  end
else
  const σ_ρ_cases_non_centered_mean = -2.2
  const σ_ρ_cases_non_centered_sd = 0.2
  const ρ_cases_init_non_centered_mean = -2.5
  const ρ_cases_init_non_centered_sd = 0.1
end

## R0
if half_R0_0
  const R₀_init_non_centered_mean = 0.0 - log(2)
else
  const R₀_init_non_centered_mean = 0.0
end
const R₀_init_non_centered_sd = 0.25
const σ_R0_non_centered_mean = -1.9
const σ_R0_non_centered_sd = 0.3

## IFR
if double_IFR_0
  const IFR_init_non_centered_mean = -5.3 + log(2)
else
  const IFR_init_non_centered_mean = -5.3
end
const IFR_init_non_centered_sd = 0.2
const σ_IFR_non_centered_mean = -2.4
const σ_IFR_non_centered_sd = 0.12

## Intitial Conditions
if half_S_0
  const S_SEI_non_centered_mean = 6 - log(2)
else
  const S_SEI_non_centered_mean = 6
end
const S_SEI_non_centered_sd = 0.5
const I_EI_non_centered_mean = 0.6
const I_EI_non_centered_sd = 0.03

## Dwell Times
# dur_latent_non_centered_mean = log(3 / 7)
# dur_latent_non_centered_sd = 0.25
# https://academic.oup.com/cid/advance-article/doi/10.1093/cid/ciab746/6359063
const dur_latent_non_centered_mean = -0.25
const dur_latent_non_centered_sd = 0.1

# dur_infectious_non_centered_mean = log(5 / 7)
# dur_infectious_non_centered_sd = 0.2
# https://bmjopen.bmj.com/content/10/8/e039856.abstract
const dur_infectious_non_centered_mean = 0.15
const dur_infectious_non_centered_sd = 0.1
const ρ_death_non_centered_mean = 2.3
const ρ_death_non_centered_sd = 0.2

## Overdispersion
include(projectdir("src/prior_constants_overdispersion.jl"))
