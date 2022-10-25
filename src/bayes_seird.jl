prob_skeleton = ODEProblem(seirdc_log_ode!,
  zeros(6),
  (0.0, obstimes[end]),
  ones(4))

@model function bayes_seird(data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, extra_ode_precision)
  # Calculate number of observed datapoints timepoints
  l_incidence = length(data_new_deaths)
  if use_seroprev
    l_prevalence = length(data_seroprev_cases)
  end

  l_param_change_times = length(param_change_times)

  if constant_R0
    l_param_change_times_R0 = 0
  else
    l_param_change_times_R0 = l_param_change_times
  end

  if constant_alpha
    l_param_change_times_alpha = 0
  else
    l_param_change_times_alpha = l_param_change_times
  end

  if constant_IFR
    l_param_change_times_IFR = 0
  else
    l_param_change_times_IFR = l_param_change_times
  end

  # Priors
  R0_params_non_centered ~ MvNormal(l_param_change_times_R0 + 2, 1) # +2, 1 for var, 1 for init
  IFR_t_params_non_centered ~ MvNormal(l_param_change_times_IFR + 2, 1) # +2, 1 for var, 1 for init

  if use_tests
    α_t_params_non_centered ~ MvNormal(l_param_change_times_alpha + 2, 1) # +2, 1 for var, 1 for init
  else
    ρ_cases_t_params_non_centered ~ MvNormal(l_param_change_times_alpha + 2, 1) # +2, 1 for var, 1 for init
  end

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

  if use_tests
    ϕ_cases_bb = exp(ϕ_cases_non_centered * ϕ_cases_bb_non_centered_sd + ϕ_cases_bb_non_centered_mean)
  else
    ϕ_cases_nb = exp(ϕ_cases_non_centered * ϕ_cases_nb_non_centered_sd + ϕ_cases_nb_non_centered_mean)
  end

  ϕ_deaths = exp(ϕ_deaths_non_centered * ϕ_deaths_non_centered_sd + ϕ_deaths_non_centered_mean)

  if !constant_R0
    σ_R0_non_centered = R0_params_non_centered[2]
    σ_R0 = exp(σ_R0_non_centered * σ_R0_non_centered_sd + σ_R0_non_centered_mean)
    log_R0_steps_non_centered = R0_params_non_centered[3:end]
  end

  if !constant_IFR
    σ_IFR_non_centered = IFR_t_params_non_centered[2]
    σ_IFR = exp(σ_IFR_non_centered * σ_IFR_non_centered_sd + σ_IFR_non_centered_mean)
    logit_IFR_t_steps_non_centered = IFR_t_params_non_centered[3:end]
  end

  if use_tests
    α_init_non_centered = α_t_params_non_centered[1]
    α_init = exp(α_init_non_centered * α_init_non_centered_sd + α_init_non_centered_mean)
  else
    ρ_cases_init_non_centered = ρ_cases_t_params_non_centered[1]
    ρ_cases_init = logistic(ρ_cases_init_non_centered * ρ_cases_init_non_centered_sd + ρ_cases_init_non_centered_mean)
  end

  if !constant_alpha
    if use_tests
      σ_α_non_centered = α_t_params_non_centered[2]
      σ_α = exp(σ_α_non_centered * σ_α_non_centered_sd + σ_α_non_centered_mean)
      log_α_t_steps_non_centered = α_t_params_non_centered[3:end]
      α_t_values_with_init = vcat(α_init, exp.(log(α_init) .+ cumsum(vec(log_α_t_steps_non_centered) * σ_α)))
    else
      σ_ρ_cases_non_centered = ρ_cases_t_params_non_centered[2]
      σ_ρ_cases = exp(σ_ρ_cases_non_centered * σ_ρ_cases_non_centered_sd + σ_ρ_cases_non_centered_mean)
      logit_ρ_cases_t_steps_non_centered = ρ_cases_t_params_non_centered[3:end]
      ρ_cases_t_values_with_init = vcat(ρ_cases_init, exp.(logit(ρ_cases_init) .+ cumsum(vec(logit_ρ_cases_t_steps_non_centered) * σ_ρ_cases)))
    end
  end

  R₀_init_non_centered = R0_params_non_centered[1]
  R₀_init = exp(R₀_init_non_centered * R₀_init_non_centered_sd + R₀_init_non_centered_mean)

  IFR_init_non_centered = IFR_t_params_non_centered[1]
  IFR_init = logistic(IFR_init_non_centered * IFR_init_non_centered_sd + IFR_init_non_centered_mean)

  β_init = R₀_init * ν

  S_SEI = logistic(S_SEI_non_centered * S_SEI_non_centered_sd + S_SEI_non_centered_mean)
  I_EI = logistic(I_EI_non_centered * I_EI_non_centered_sd + I_EI_non_centered_mean)
  S_init = S_SEI * popsize
  I_init = max(I_EI * (popsize - S_init), 1) # Make sure at least 1 Infectious
  E_init = max(popsize - (S_init + I_init), 1) # Make sure at least 1 Exposed
  u0 = [S_init, E_init, I_init, 1.0, 1.0, I_init] # Intialize with 1 in R and D so there are no problems when we log for ODE

  # Time-varying parameters
  if !constant_R0
    β_t_values_no_init = exp.(log(R₀_init) .+ cumsum(vec(log_R0_steps_non_centered) * σ_R0)) * ν
  end

  if !constant_IFR
    IFR_t_values_no_init = logistic.(logit(IFR_init) .+ cumsum(vec(logit_IFR_t_steps_non_centered) * σ_IFR))
  end

  function param_affect_β_IFR!(integrator)
    ind_t = searchsortedfirst(param_change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
    if !constant_R0
      integrator.p[1] = β_t_values_no_init[ind_t] # Replace β with a new value from β_t_values
    end
    if !constant_IFR
      integrator.p[4] = IFR_t_values_no_init[ind_t] # Replace IFR with a new value from IFR_t_values
    end
  end
  param_callback = PresetTimeCallback(param_change_times, param_affect_β_IFR!, save_positions = (false, false))

  prob = remake(prob_skeleton,
    u0 = log.(u0),
    p = [β_init, γ, ν, IFR_init],
    tspan = (0.0, obstimes[end]))

  if extra_ode_precision
    sol = solve(prob, Tsit5(), callback = param_callback, saveat = obstimes, save_start = true, verbose = false, abstol = 1e-11, reltol = 1e-8)
  else
    sol = solve(prob, Tsit5(), callback = param_callback, saveat = obstimes, save_start = true, verbose = false, abstol = 1e-9, reltol = 1e-6)
  end

  # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
  if sol.retcode != :Success
    Turing.@addlogprob! -Inf
    return
  end

  sol_reg_scale_array = exp.(Array(sol))

  sol_new_deaths = sol_reg_scale_array[5, 2:end] - sol_reg_scale_array[5, 1:(end-1)]
  sol_new_cases = sol_reg_scale_array[6, 2:end] - sol_reg_scale_array[6, 1:(end-1)]

  deaths_mean = ρ_death * sol_new_deaths

  if constant_alpha
    if use_tests
      cases_bb_mean = logistic.(α_init .+ logit.(sol_new_cases / popsize))
    else
      cases_nb_mean = ρ_cases_init .* sol_new_cases
    end
  else
    if use_tests
      cases_bb_mean = logistic.(α_t_values_with_init .+ logit.(sol_new_cases / popsize))
    else
      cases_nb_mean = ρ_cases_t_values_with_init .* sol_new_cases
    end
  end

  if use_seroprev
    compartments_at_seroprev_times = [exp.(sol.u[findfirst(sol.t .== seroprev_time)]) for seroprev_time in seroprev_times] # Collect S, E, I, R, D, C at seroprev_times
    seroprev_mean = [compartments[4] / sum(compartments[1:4]) for compartments in compartments_at_seroprev_times]
  end

  for i in 1:l_incidence
    data_new_deaths[i] ~ NegativeBinomial2(max(deaths_mean[i], 0.0), ϕ_deaths) # In exploratory phase of MCMC, sometimes you get weird numerical errors
    if use_tests
      data_new_cases[i] ~ BetaBinomial2(tests[i], max(cases_bb_mean[i], 0.0), ϕ_cases_bb)
    else
      data_new_cases[i] ~ NegativeBinomial2(max(cases_nb_mean[i], 0.0), ϕ_cases_nb)
    end
  end

  if use_seroprev
    for i in 1:l_prevalence
      data_seroprev_cases[i] ~ Binomial(seroprev_tests[i], seroprev_mean[i])
    end
  end

  # Generated quantities
  gq = (
    S_SEI = S_SEI,
    I_EI = I_EI,
    S = sol_reg_scale_array[1, :],
    E = sol_reg_scale_array[2, :],
    I = sol_reg_scale_array[3, :],
    R = sol_reg_scale_array[4, :],
    D = sol_reg_scale_array[5, :],
    C = sol_reg_scale_array[6, :],
    prevalence = sol_reg_scale_array[2, :] + sol_reg_scale_array[3, :],
    dur_latent_days = 7 / γ,
    dur_infectious_days = 7 / ν,
    ϕ_deaths = ϕ_deaths,
    ρ_death = ρ_death,
    deaths_mean = deaths_mean)
  if constant_alpha
    if use_tests
      gq = merge(gq, (α = α_init,))
    else
      gq = merge(gq, (ρ_cases = ρ_cases_init,))
    end
  else
    if use_tests
      gq = merge(gq, (σ_α = σ_α, α_t = α_t_values_with_init, ϕ_cases_bb = ϕ_cases_bb, cases_bb_mean = cases_bb_mean, cases_mean = tests .* cases_bb_mean))
    else
      gq = merge(gq, (σ_ρ_cases = σ_ρ_cases, ρ_cases_t = ρ_cases_t_values_with_init, ϕ_cases_nb = ϕ_cases_nb, cases_nb_mean = cases_nb_mean, cases_mean = cases_nb_mean))
    end
  end

  if constant_R0
    β = β_init
    R₀ = β / ν
    Rₜ_t = R₀ * gq[:S][1:(end-1)] / popsize
    gq = merge(gq, (β = β, R₀ = R₀, Rₜ_t = Rₜ_t))
  else
    β_t = vcat(β_init, β_t_values_no_init)
    R₀_t = β_t / ν
    Rₜ_t = R₀_t .* gq[:S][1:(end-1)] / popsize
    gq = merge(gq, (β_t = β_t, R₀_t = R₀_t, Rₜ_t = Rₜ_t, σ_R0 = σ_R0))
  end

  if constant_IFR
    gq = merge(gq, (IFR = IFR_init,))
  else
    gq = merge(gq, (IFR_t = vcat(IFR_init, IFR_t_values_no_init), σ_IFR = σ_IFR))
  end

  if use_seroprev
    gq = merge(gq, (seroprev_mean = seroprev_mean,))
  end

  return gq
end