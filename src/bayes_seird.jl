prob = ODEProblem{true}(seirdc_log_ode!,
  zeros(6),
  (0.0, obstimes[end]),
  ones(4))

@model function bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests::Bool, use_seroprev::Bool, constant_R0::Bool, constant_alpha::Bool, constant_IFR::Bool, DEAlgorithm::DEAlgorithm, abstol, reltol)
  l_incidence = length(data_new_deaths)

  if use_seroprev
    l_prevalence = length(data_seroprev_cases)
  end

  l_param_change_times = length(param_change_times)
  l_param_change_times_R0 = constant_R0 ? 0 : l_param_change_times
  l_param_change_times_alpha = constant_alpha ? 0 : l_param_change_times
  l_param_change_times_IFR = constant_IFR ? 0 : l_param_change_times

  # Priors
  R0_params_non_centered ~ MvNormal(Zeros(l_param_change_times_R0 + 2), I) # +2, 1 for var, 1 for init
  IFR_t_params_non_centered ~ MvNormal(Zeros(l_param_change_times_IFR + 2), I) # +2, 1 for var, 1 for init

  if use_tests
    α_t_params_non_centered ~ MvNormal(Zeros(l_param_change_times_alpha + 2), I) # +2, 1 for var, 1 for init
  else
    ρ_cases_t_params_non_centered ~ MvNormal(Zeros(l_param_change_times_alpha + 2), I) # +2, 1 for var, 1 for init
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

  R₀_init_non_centered = R0_params_non_centered[1]
  R₀_init = exp(R₀_init_non_centered * R₀_init_non_centered_sd + R₀_init_non_centered_mean)
  β_init = R₀_init * ν
  β_t_values_no_init = constant_R0 ? repeat([β_init], l_param_change_times) : exp.(log(R₀_init) .+ cumsum(log_R0_steps_non_centered) * σ_R0) * ν

  if !constant_IFR
    σ_IFR_non_centered = IFR_t_params_non_centered[2]
    σ_IFR = exp(σ_IFR_non_centered * σ_IFR_non_centered_sd + σ_IFR_non_centered_mean)
    logit_IFR_t_steps_non_centered = IFR_t_params_non_centered[3:end]
  end

  IFR_init_non_centered = IFR_t_params_non_centered[1]
  IFR_init = logistic(IFR_init_non_centered * IFR_init_non_centered_sd + IFR_init_non_centered_mean)
  IFR_t_values_no_init = constant_IFR ? repeat([IFR_init], l_param_change_times) : logistic.(logit(IFR_init) .+ cumsum(logit_IFR_t_steps_non_centered) * σ_IFR)

  if !constant_alpha
    if use_tests
      σ_α_non_centered = α_t_params_non_centered[2]
      σ_α = exp(σ_α_non_centered * σ_α_non_centered_sd + σ_α_non_centered_mean)
      log_α_t_steps_non_centered = α_t_params_non_centered[3:end]
    else
      σ_ρ_cases_non_centered = ρ_cases_t_params_non_centered[2]
      σ_ρ_cases = exp(σ_ρ_cases_non_centered * σ_ρ_cases_non_centered_sd + σ_ρ_cases_non_centered_mean)
      logit_ρ_cases_t_steps_non_centered = ρ_cases_t_params_non_centered[3:end]
    end
  end

  if use_tests
    α_init_non_centered = α_t_params_non_centered[1]
    α_init = exp(α_init_non_centered * α_init_non_centered_sd + α_init_non_centered_mean)
    α_t_values_with_init = constant_alpha ? repeat([α_init], l_param_change_times + 1) : vcat(α_init, exp.(log(α_init) .+ cumsum(vec(log_α_t_steps_non_centered) * σ_α)))
  else
    ρ_cases_init_non_centered = ρ_cases_t_params_non_centered[1]
    ρ_cases_init = logistic(ρ_cases_init_non_centered * ρ_cases_init_non_centered_sd + ρ_cases_init_non_centered_mean)
    ρ_cases_t_values_with_init = constant_alpha ? repeat([ρ_cases_init], l_param_change_times + 1) : vcat(ρ_cases_init, exp.(logit(ρ_cases_init) .+ cumsum(vec(logit_ρ_cases_t_steps_non_centered) * σ_ρ_cases)))
  end

  # Initial Conditions
  S_SEI = logistic(S_SEI_non_centered * S_SEI_non_centered_sd + S_SEI_non_centered_mean)
  I_EI = logistic(I_EI_non_centered * I_EI_non_centered_sd + I_EI_non_centered_mean)
  S_init = S_SEI * popsize
  I_init = max(I_EI * (popsize - S_init), 1) # Make sure at least 1 Infectious
  E_init = max(popsize - (S_init + I_init), 1) # Make sure at least 1 Exposed
  u0 = [S_init, E_init, I_init, 1.0, 1.0, I_init] # Intialize with 1 in R and D so there are no problems when we log for ODE

  function param_affect_β_IFR!(integrator)
    ind_t = searchsortedfirst(param_change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
    integrator.p[1] = β_t_values_no_init[ind_t] # Replace β with a new value from β_t_values
    integrator.p[4] = IFR_t_values_no_init[ind_t] # Replace IFR with a new value from IFR_t_values
  end
  param_callback = PresetTimeCallback(param_change_times, param_affect_β_IFR!, save_positions=(false, false))

  sol = solve(prob, DEAlgorithm; callback=param_callback, saveat=obstimes, save_start=true, verbose=false, abstol=abstol, reltol=reltol, u0=log.(u0), p=[β_init, γ, ν, IFR_init], tspan=(0.0, obstimes[end]))

  # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
  if sol.retcode != :Success
    Turing.@addlogprob! -Inf
    return
  end

  sol_reg_scale_array = exp.(Array(sol))

  sol_new_deaths = sol_reg_scale_array[5, 2:end] - sol_reg_scale_array[5, 1:(end-1)]
  sol_new_cases = sol_reg_scale_array[6, 2:end] - sol_reg_scale_array[6, 1:(end-1)]

  deaths_mean = max.(ρ_death * sol_new_deaths, 0.0)

  if use_tests
    cases_bb_mean = max.(logistic.(α_t_values_with_init .+ logit.(sol_new_cases / popsize)), 0.0)
  else
    cases_nb_mean = max.(ρ_cases_t_values_with_init .* sol_new_cases, 0.0)
  end

  if use_seroprev
    compartments_at_seroprev_times = sol_reg_scale_array[:, findall(x -> x ∈ seroprev_times, sol.t)] # Collect S, E, I, R, D, C at seroprev_times
    seroprev_mean = compartments_at_seroprev_times[4, :] ./ vec(sum(compartments_at_seroprev_times[1:4, :], dims=1))
  end

  for i in 1:l_incidence
    data_new_deaths[i] ~ NegativeBinomial2(deaths_mean[i], ϕ_deaths) # In exploratory phase of MCMC, sometimes you get weird numerical errors
    if use_tests
      data_new_cases[i] ~ BetaBinomial2(tests[i], cases_bb_mean[i], ϕ_cases_bb)
    else
      data_new_cases[i] ~ NegativeBinomial2(cases_nb_mean[i], ϕ_cases_nb)
    end
  end

  if use_seroprev
    for i in 1:l_prevalence
      data_seroprev_cases[i] ~ Binomial(seroprev_tests[i], seroprev_mean[i])
    end
  end

  S = sol_reg_scale_array[1, :]
  β_t = vcat(β_init, β_t_values_no_init)
  R₀_t = β_t / ν

  # Generated quantities
  gq_base = (
    S_SEI=S_SEI,
    I_EI=I_EI,
    S=S,
    E=sol_reg_scale_array[2, :],
    I=sol_reg_scale_array[3, :],
    R=sol_reg_scale_array[4, :],
    D=sol_reg_scale_array[5, :],
    C=sol_reg_scale_array[6, :],
    prevalence=sol_reg_scale_array[2, :] + sol_reg_scale_array[3, :],
    dur_latent_days=time_interval_in_days / γ,
    dur_infectious_days=time_interval_in_days / ν,
    ϕ_deaths=ϕ_deaths,
    ρ_death=ρ_death,
    deaths_mean=deaths_mean,
    β_t=β_t,
    R₀_t=R₀_t,
    Rₜ_t=R₀_t .* S[1:(end-1)] / popsize,
    IFR_t=vcat(IFR_init, IFR_t_values_no_init))

  if use_tests
    gq = merge(gq_base, (α_t=α_t_values_with_init, ϕ_cases_bb=ϕ_cases_bb, cases_bb_mean=cases_bb_mean, cases_mean=tests .* cases_bb_mean))
  else
    gq = merge(gq_base, (ρ_cases_t=ρ_cases_t_values_with_init, ϕ_cases_nb=ϕ_cases_nb, cases_nb_mean=cases_nb_mean, cases_mean=cases_nb_mean))
  end

  if !constant_R0
    gq = merge(gq, (σ_R0=σ_R0,))
  end

  if !constant_IFR
    gq = merge(gq, (σ_IFR=σ_IFR,))
  end

  if !constant_alpha
    if use_tests
      gq = merge(gq, (σ_α=σ_α,))
    else
      gq = merge(gq, (σ_ρ_cases=σ_ρ_cases,))
    end
  end

  if use_seroprev
    gq = merge(gq, (seroprev_mean=seroprev_mean,))
  end

  return gq
end