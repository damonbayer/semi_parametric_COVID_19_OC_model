using DifferentialEquations
using LinearAlgebra
using Turing
using FillArrays
using Random
using LogExpFunctions
using JLD2
using FileIO
using DrWatson
using CSV
using DataFrames
using semi_parametric_COVID_19_OC_model
max_t = 20.0
obstimes = collect(1:max_t)
function seirdc_age_log_ode(du, u, p, t)
    (S_o, E_o, I_o, R_o, D_o, C_o, S_y, E_y, I_y, R_y, D_y, C_y) = exp.(u)
    (β_oo, γ, ν, IFR_o, β_yy, IFR_y, β_yo, β_oy) = p

    N = S_o + E_o + I_o + R_o + D_o + S_y + E_y + I_y + R_y + D_y

    infection_oo = β_oo * I_o * S_o / N
    infection_oy = β_oy * I_o * S_y / N
    infection_yo = β_yo * I_y * S_o / N
    infection_yy = β_yy * I_y * S_y / N

    infection_o = infection_oo + infection_yo
    infection_y = infection_yy + infection_oy

    progression_o = γ * E_o
    progression_y = γ * E_y

    recovery_o = ν * (1 - IFR_o) * I_o
    recovery_y = ν * (1 - IFR_y) * I_y

    death_o = ν * IFR_o * I_o
    death_y = ν * IFR_y * I_y

    @inbounds begin
        du[1] = -infection_o / S_o # S
        du[2] = (infection_o - progression_o) / E_o # E
        du[3] = (progression_o - recovery_o - death_o) / I_o # I
        du[4] = recovery_o / R_o # R
        du[5] = death_o / D_o # D
        du[6] = progression_o / C_o # Cumulative Progressions
        du[7] = -infection_y / S_y # S
        du[8] = (infection_y - progression_y) / E_y # E
        du[9] = (progression_y - recovery_y - death_y) / I_y # I
        du[10] = recovery_y / R_y # R
        du[11] = death_y / D_y # D
        du[12] = progression_y / C_y # Cumulative Progressions
    end
    nothing
end

popsize_y = 10_000
E_y_init = 40
I_y_init = 20
R_y_init = 1
D_y_init = 1
C_y_init = I_y_init + R_y_init + D_y_init
S_y_init = popsize_y - sum([E_y_init, I_y_init, R_y_init, D_y_init])

popsize_o = 1_000
E_o_init = 40
I_o_init = 20
R_o_init = 1
D_o_init = 1
C_o_init = I_o_init + R_o_init + D_o_init
S_o_init = popsize_o - sum([E_o_init, I_o_init, R_o_init, D_o_init])

const S_init = S_y_init + S_o_init
const I_init = I_y_init + I_o_init
const E_init = E_y_init + E_o_init
const popsize = popsize_y + popsize_o

ρ = popsize_o / popsize

init_state = [S_o_init, E_o_init, I_o_init, R_o_init, D_o_init, C_o_init,
    S_y_init, E_y_init, I_y_init, R_y_init, D_y_init, C_y_init]

β_oo = 4
γ = 1
ν = 2
IFR_o = 0.1
β_yy = β_oo
IFR_y = 0.01
β_yo = β_oy = 1.0

NGM = [(ρ * β_oo) / ν (ρ * β_oy) / ν; ((1 - ρ) * β_yo) / ν ((1 - ρ) * β_yy) / ν]

prob = ODEProblem(seirdc_age_log_ode,
    log.(init_state),
    (0.0, max_t),
    (β_oo, γ, ν, IFR_o, β_yy, IFR_y, β_yo, β_oy))

sol = solve(prob, Tsit5(), saveat=1.0, save_start=true, verbose=false, abstol=1e-11, reltol=1e-8);
sol_array = exp.(Array(sol))

new_latent_deaths_o = sol_array[5, 2:end] - sol_array[5, 1:(end-1)]
new_latent_deaths_y = sol_array[11, 2:end] - sol_array[11, 1:(end-1)]
new_latent_deaths_total = new_latent_deaths_o + new_latent_deaths_y

new_latent_cases_o = sol_array[6, 2:end] - sol_array[6, 1:(end-1)]
new_latent_cases_y = sol_array[12, 2:end] - sol_array[12, 1:(end-1)]
new_latent_cases_total = new_latent_cases_o + new_latent_cases_y

Random.seed!(2)
data_new_cases = vcat(rand.(Poisson.(new_latent_cases_total), 1)...)
Random.seed!(2)
data_new_deaths = vcat(rand.(Poisson.(new_latent_deaths_total), 1)...)

dat = DataFrame(time = Int.(obstimes), cases = data_new_cases, deaths = data_new_deaths)

CSV.write(projectdir("illustrative_examples", "age_structure", "data", "data.csv"), dat)

true_generated_quantities = rename!(DataFrame(transpose(sol_array), :auto), [:S_o, :E_o, :I_o, :R_o, :D_o, :C_o, :S_y, :E_y, :I_y, :R_y, :D_y, :C_y])
true_generated_quantities.time = vcat(0, obstimes)
true_generated_quantities.S_total = true_generated_quantities[!, :S_o] + true_generated_quantities[!, :S_y]
true_generated_quantities.E_total = true_generated_quantities[!, :E_o] + true_generated_quantities[!, :E_y]
true_generated_quantities.I_total = true_generated_quantities[!, :I_o] + true_generated_quantities[!, :I_y]
true_generated_quantities.R_total = true_generated_quantities[!, :R_o] + true_generated_quantities[!, :R_y]
true_generated_quantities.D_total = true_generated_quantities[!, :D_o] + true_generated_quantities[!, :D_y]

true_generated_quantities.latent_deaths_o = vcat(missing, new_latent_deaths_o)
true_generated_quantities.latent_deaths_y = vcat(missing, new_latent_deaths_y)
true_generated_quantities.latent_deaths_total = vcat(missing, new_latent_deaths_total)
true_generated_quantities.latent_cases_o = vcat(missing, new_latent_cases_o)
true_generated_quantities.latent_cases_y = vcat(missing, new_latent_cases_y)
true_generated_quantities.latent_cases_total = vcat(missing, new_latent_cases_total)

true_generated_quantities.dur_latent_days=vcat(7 / γ, repeat([missing], nrow(true_generated_quantities)-1))
true_generated_quantities.dur_infectious_days=vcat(7 / ν, repeat([missing], nrow(true_generated_quantities)-1))
true_generated_quantities.deaths_mean=vcat(missing, new_latent_deaths_o + new_latent_deaths_y)
true_generated_quantities.cases_mean=vcat(missing, new_latent_cases_o + new_latent_cases_y)
true_generated_quantities.R₀=vcat(maximum(eigvals(NGM)), repeat([missing], nrow(true_generated_quantities)-1))
true_generated_quantities.β = ν * true_generated_quantities.R₀
true_generated_quantities.Rₜ_t = true_generated_quantities.R₀[1] .* true_generated_quantities.S_total / popsize
true_generated_quantities.IFR_t = (true_generated_quantities[!, :I_o] * IFR_o + true_generated_quantities[!, :I_y] * IFR_y) ./ (true_generated_quantities[!, :I_o] + true_generated_quantities[!, :I_y])

CSV.write(projectdir("illustrative_examples", "age_structure", "data", "true_generated_quantities.csv"), true_generated_quantities)

function seirdc_log_ode!(du, u, p, t)
    S, E, I, R, D, C = exp.(u)
    β, γ, ν, IFR = p
    N = S + E + I + R + D

    infection = β * I * S / N
    progression = γ * E
    recovery = ν * (1 - IFR) * I
    death = ν * IFR * I

    @inbounds begin
        du[1] = -infection / S # S
        du[2] = (infection - progression) / E # E
        du[3] = (progression - (recovery + death)) / I # I
        du[4] = recovery / R # R
        du[5] = death / D # D
        du[6] = progression / C # Cumulative Progressions
    end
    nothing
end


# Can you infer it?
obstimes = collect(1:max_t)
param_change_times = obstimes[1:(end-1)]

prob = ODEProblem{true}(seirdc_log_ode!,
    zeros(6),
    (0.0, obstimes[end]),
    ones(4))

@model function bayes_seird(prob, data_new_deaths, data_new_cases, obstimes, param_change_times, extra_ode_precision)
    # Calculate number of observed datapoints timepoints
    l_incidence = length(data_new_deaths)
    l_param_change_times = length(param_change_times)

    # Priors
    R₀_init_non_centered ~ Normal()
    IFR_t_params_non_centered ~ MvNormal(Zeros(l_param_change_times + 2), I) # +2, 1 for var, 1 for init
    dur_latent_non_centered ~ Normal()
    dur_infectious_non_centered ~ Normal()

    # Transformations
    γ = exp(-(dur_latent_non_centered * 0.2 + 0.25))

    ν = exp(-(dur_infectious_non_centered * 0.1 - 1))

    σ_IFR_non_centered = IFR_t_params_non_centered[2]
    σ_IFR = exp(σ_IFR_non_centered * 0.13 + -1.5)
    logit_IFR_t_steps_non_centered = IFR_t_params_non_centered[3:end]

    R₀_init = exp(R₀_init_non_centered * 0.1 + 0.5)

    IFR_init_non_centered = IFR_t_params_non_centered[1]

    IFR_init = logistic(IFR_init_non_centered * 0.3 - 2.6)


    β_init = R₀_init * ν

    u0 = [S_init, E_init, I_init, 1.0, 1.0, E_init + I_init] # Intialize with 1 in R and D so there are no problems when we log for ODE
    IFR_t_values_no_init = logistic.(logit(IFR_init) .+ cumsum(logit_IFR_t_steps_non_centered) * σ_IFR)

    function param_affect_β_IFR!(integrator)
        ind_t = searchsortedfirst(param_change_times, integrator.t) # Find the index of param_change_times that contains the current timestep
        integrator.p[4] = IFR_t_values_no_init[ind_t] # Replace IFR with a new value from IFR_t_values
    end

    param_callback = PresetTimeCallback(param_change_times, param_affect_β_IFR!, save_positions=(false, false))

    if extra_ode_precision
        sol = solve(prob, Tsit5(); callback=param_callback, saveat=obstimes, save_start=true, verbose=false, abstol=1e-10, reltol=1e-7, u0=log.(u0), p=[β_init, γ, ν, IFR_init])
    else
        sol = solve(prob, Tsit5(); callback=param_callback, saveat=obstimes, save_start=true, verbose=false, abstol=1e-8, reltol=1e-5, u0=log.(u0), p=[β_init, γ, ν, IFR_init])
    end

    # If the ODE solver fails, reject the sample by adding -Inf to the likelihood
    if sol.retcode != :Success
        Turing.@addlogprob! -Inf
        return
    end

    sol_reg_scale_array = exp.(Array(sol))

    sol_new_deaths = sol_reg_scale_array[5, 2:end] - sol_reg_scale_array[5, 1:(end-1)]
    sol_new_cases = sol_reg_scale_array[6, 2:end] - sol_reg_scale_array[6, 1:(end-1)]

    deaths_mean = sol_new_deaths
    cases_mean = sol_new_cases

    for i in 1:l_incidence
        data_new_deaths[i] ~ Poisson(max(deaths_mean[i], 0.0)) # In exploratory phase of MCMC, sometimes you get weird numerical errors
        data_new_cases[i] ~ Poisson(max(cases_mean[i], 0.0))
    end

    # # Generated quantities
    gq = (
        S=sol_reg_scale_array[1, :],
        E=sol_reg_scale_array[2, :],
        I=sol_reg_scale_array[3, :],
        R=sol_reg_scale_array[4, :],
        D=sol_reg_scale_array[5, :],
        dur_latent_days=7 / γ,
        dur_infectious_days=7 / ν,
        deaths_mean=deaths_mean,
        cases_mean=cases_mean,
        β=β_init,
        R₀=β_init / ν,
        Rₜ_t=β_init / ν * sol_reg_scale_array[1, :][1:(end-1)] / popsize,
        IFR_t=vcat(IFR_init, IFR_t_values_no_init),
        σ_IFR=σ_IFR)
    return gq
end

my_model = bayes_seird(prob, data_new_deaths, data_new_cases, obstimes, param_change_times, false)
my_model_missing = bayes_seird(prob, repeat([missing], length(data_new_deaths)), repeat([missing], length(data_new_cases)), obstimes, param_change_times, true)

mkpath(projectdir("illustrative_examples", "age_structure"))

# Sample Prior
Random.seed!(1)
prior_samples = sample(my_model, Prior(), 1000)
wsave(projectdir("illustrative_examples", "age_structure", "data", "prior_samples.jld2"), @dict prior_samples)

Random.seed!(1)
prior_predictive = predict(my_model_missing, prior_samples)
CSV.write(projectdir("illustrative_examples", "age_structure", "data", "prior_predictive.csv"), DataFrame(prior_predictive))

Random.seed!(1)
prior_generated_quantities = hcat(Chains(generated_quantities(my_model, prior_samples)), getlogp_chain(prior_samples))
CSV.write(projectdir("illustrative_examples", "age_structure", "data", "prior_generated_quantities.csv"), DataFrame(prior_generated_quantities))

# Sample Posterior
Random.seed!(1)
posterior_samples = sample(my_model, NUTS(), MCMCThreads(), 500, 4)
wsave(projectdir("illustrative_examples", "age_structure", "data", "posterior_samples.jld2"), @dict posterior_samples)

Random.seed!(1)
posterior_predictive = predict(my_model_missing, posterior_samples)
CSV.write(projectdir("illustrative_examples", "age_structure", "data", "posterior_predictive.csv"), DataFrame(posterior_predictive))

Random.seed!(1)
posterior_generated_quantities = hcat(Chains(generated_quantities(my_model, posterior_samples)), getlogp_chain(posterior_samples))
CSV.write(projectdir("illustrative_examples", "age_structure", "data", "posterior_generated_quantities.csv"), DataFrame(posterior_generated_quantities))