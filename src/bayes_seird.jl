@model function bayes_seird(data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times)
  # Calculate number of observed datapoints timepoints
  l_incidence = length(data_new_deaths)
  l_prevalence = length(data_seroprev_cases)
  l_param_change_times = length(param_change_times)

  # Priors
  R0_params_non_centered ~ MvNormal(l_param_change_times + 2, 1) # +2, 1 for var, 1 for init
  IFR_t_params_non_centered ~ MvNormal(l_param_change_times + 2, 1) # +2, 1 for var, 1 for init
  α_t_params_non_centered ~ MvNormal(l_param_change_times + 2, 1) # +2, 1 for var, 1 for init
  S_SEI_non_centered ~ Normal()
  I_EI_non_centered ~ Normal()
  dur_latent_non_centered ~ Normal()
  dur_infectious_non_centered ~ Normal()
  ϕ_cases_non_centered ~ Normal()
  ϕ_deaths_non_centered ~ Normal()
  ρ_death_non_centered ~ Normal()

  # Transformations
  γ = exp(-(dur_latent_non_centered * dur_latent_non_centered_sd + dur_latent_non_centered_mean))
  ν = exp(-(dur_infectious_non_centered * dur_infectious_non_centered_sd + dur_infectious_non_centered_mean))
  ρ_death = logistic(ρ_death_non_centered * ρ_death_non_centered_sd + ρ_death_non_centered_mean)

  ϕ_cases = exp(ϕ_cases_non_centered * 0.24 + 7.0)
  ϕ_deaths = exp(ϕ_deaths_non_centered * 0.46 + 5.4)

  σ_R0_non_centered = R0_params_non_centered[1]
  σ_IFR_non_centered = IFR_t_params_non_centered[1]
  σ_α_non_centered = α_t_params_non_centered[1]

  σ_R0 = exp(σ_R0_non_centered * σ_R0_non_centered_sd + σ_R0_non_centered_mean)
  σ_IFR = exp(σ_IFR_non_centered * σ_IFR_non_centered_sd + σ_IFR_non_centered_mean)
  σ_α = exp(σ_α_non_centered * σ_α_non_centered_sd + σ_α_non_centered_mean)

  R₀_init_non_centered = R0_params_non_centered[2]
  IFR_init_non_centered = IFR_t_params_non_centered[2]
  α_init_non_centered = α_t_params_non_centered[2]

  R₀_init = exp(R₀_init_non_centered * R₀_init_non_centered_sd + R₀_init_non_centered_mean)
  β_init = R₀_init * ν
  IFR_init = logistic(IFR_init_non_centered * IFR_init_non_centered_sd + IFR_init_non_centered_mean)
  α_init = exp(α_init_non_centered * α_init_non_centered_sd + α_init_non_centered_mean)

  log_R0_steps_non_centered = R0_params_non_centered[3:end]
  logit_IFR_t_steps_non_centered = IFR_t_params_non_centered[3:end]
  log_α_t_steps_non_centered = α_t_params_non_centered[3:end]

  S_SEI = logistic(S_SEI_non_centered * S_SEI_non_centered_sd + S_SEI_non_centered_mean)
  I_EI = logistic(I_EI_non_centered * I_EI_non_centered_sd + I_EI_non_centered_mean)

  S_init = S_SEI * popsize
  I_init = max(I_EI * (popsize - S_init), 1) # Make sure at least 1 Infectious
  E_init = max(popsize - (S_init + I_init), 1) # Make sure at least 1 Exposed
  u0 = [S_init, E_init, I_init, 1.0, 1.0, I_init] # Intialize with 1 in R and D so there are no problems when we log for ODE

  # Time-varying parameters
  β_t_values_no_init = exp.(log(R₀_init) .+ cumsum(vec(log_R0_steps_non_centered) * σ_R0)) * ν
  β_t_values_with_init = vcat(β_init, β_t_values_no_init)

  IFR_t_values_no_init = logistic.(logit(IFR_init) .+ cumsum(vec(logit_IFR_t_steps_non_centered) * σ_IFR))
  IFR_t_values_with_init = vcat(IFR_init, IFR_t_values_no_init)

  α_t_values_with_init = vcat(α_init, exp.(log(α_init) .+ cumsum(vec(log_α_t_steps_non_centered) * σ_α)))

  prob = ODEProblem(seirdc_log_ode!,
    log.(u0),
    (0.0, obstimes[end]),
    [β_init, γ, ν, IFR_init])

  function param_affect_β_IFR!(integrator)
    ind_t = searchsortedfirst(param_change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
    integrator.p[1] = β_t_values_no_init[ind_t] # Replace β with a new value from β_t_values
    integrator.p[4] = IFR_t_values_no_init[ind_t] # Replace IFR with a new value from IFR_t_values
  end

  param_callback = PresetTimeCallback(param_change_times, param_affect_β_IFR!, save_positions = (false, false))

  # Solve the ODE at intervals of 1.0, could also solve at obstimes
  sol = solve(prob, Tsit5(), callback = param_callback, saveat = obstimes, save_start = true, verbose = false)
  
  # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
  if sol.retcode != :Success
    Turing.@addlogprob! -Inf
    return
  end

  sol_reg_scale_array = exp.(Array(sol))

  sol_new_deaths = sol_reg_scale_array[5, 2:end] - sol_reg_scale_array[5, 1:(end-1)]
  sol_new_cases = sol_reg_scale_array[6, 2:end] - sol_reg_scale_array[6, 1:(end-1)]

  deaths_mean = ρ_death * sol_new_deaths
  # cases_pos_mean = (exp.(α_t_values_with_init) .* sol_new_cases / popsize) ./ (expm1.(α_t_values_with_init) .* sol_new_cases / popsize .+ 1)
  cases_pos_mean = logistic.(α_t_values_with_init .+ logit.(sol_new_cases / popsize))

  compartments_at_seroprev_times = [exp.(sol.u[findfirst(sol.t .== seroprev_time)]) for seroprev_time in seroprev_times] # Collect S, E, I, R, D, C at seroprev_times
  seroprev = [compartments[4] / sum(compartments[1:4]) for compartments in compartments_at_seroprev_times]
  
  for i in 1:l_incidence
    data_new_deaths[i] ~ NegativeBinomial2(max(deaths_mean[i], 0.0), ϕ_deaths) # In exploratory phase of MCMC, sometimes you get weird numerical errors
    data_new_cases[i] ~ BetaBinomial2(tests[i], max(cases_pos_mean[i], 0.0), ϕ_cases)
  end

  for i in 1:l_prevalence
    data_seroprev_cases[i] ~ Binomial(seroprev_tests[i], seroprev[i])
  end

  # Generated quantities
  S = sol_reg_scale_array[1, :]
  R₀_t_values_with_init = β_t_values_with_init / ν
  Rₜ_t_values_with_init = R₀_t_values_with_init .* S[1:(end-1)] / popsize

  return (
    γ = γ,
    ν = ν,
    ϕ_cases = ϕ_cases,
    ϕ_deaths = ϕ_deaths,
    ρ_death = ρ_death,
    β_t_values = β_t_values_with_init,
    R₀_t_values = R₀_t_values_with_init,
    Rₜ_t_values = Rₜ_t_values_with_init,
    IFR_t_values = IFR_t_values_with_init,
    α_t_values = α_t_values_with_init,
    S_SEI,
    I_EI,
    σ_R0,
    σ_IFR,
    σ_α,
    S = sol_reg_scale_array[1, :],
    E = sol_reg_scale_array[2, :],
    I = sol_reg_scale_array[3, :],
    R = sol_reg_scale_array[4, :],
    D = sol_reg_scale_array[5, :],
    deaths_mean = deaths_mean,
    cases_pos_mean = cases_pos_mean,
    cases_mean = tests .* cases_pos_mean,
    seroprev = seroprev
  )
end
