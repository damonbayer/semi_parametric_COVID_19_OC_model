seed = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])

using Revise
using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using LogExpFunctions
using Turing
using Random
using JLD2
using FileIO
using DrWatson
using semi_parametric_COVID_19_OC_model

max_t = 42.0
use_tests = true
use_seroprev = true
constant_R0 = false
constant_alpha = false
constant_IFR = false
double_IFR_0 = false
half_alpha_0 = false
half_S_0 = false
half_R0_0 = false

simulated_dict = @dict seed max_t use_tests use_seroprev constant_R0 constant_alpha constant_IFR double_IFR_0 half_alpha_0 half_S_0 half_R0_0

## Control Parameters
n_samples = 10_000

## Define Prior Constants
include(projectdir("src/prior_constants.jl"))

## Define Model
Turing.setadbackend(:forwarddiff)
include(projectdir("src/seirdc_log_ode.jl"))

## Load Data 
include(projectdir("src/load_process_data.jl"))
include(projectdir("src/bayes_seird.jl"))

## Create Models
my_model = bayes_seird(data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, false)
my_model_forecast_missing = bayes_seird(missing_new_deaths_forecast, missing_new_cases_forecast, tests_forecast, missing_seroprev_cases_forecast, seroprev_tests_forecast, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, true)

MAP_init = optimize_many_MAP(my_model, 10, 1, true)[1]
MAP_chain = Chains([MAP_init], flattened_varnames_list(my_model))


Random.seed!(seed)
data_chain = predict(my_model_forecast_missing, MAP_chain)


mkpath(projectdir("data", "simulated_data"))
CSV.write(projectdir("data", "simulated_data", savename("simulated_data", simulated_dict, "csv")), DataFrame(data_chain))

simulated_data_new_deaths = vcat(get(data_chain, :data_new_deaths)[:data_new_deaths]...)
simulated_data_new_cases = vcat(get(data_chain, :data_new_cases)[:data_new_cases]...)
simulated_data_seroprev_cases = vcat(get(data_chain, :data_seroprev_cases)[:data_seroprev_cases]...)

my_model_simulated = bayes_seird(simulated_data_new_deaths, simulated_data_new_cases, tests_forecast, simulated_data_seroprev_cases, seroprev_tests_forecast, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, false)

alg = Gibbs(NUTS(-1, 0.65, :dur_latent_non_centered, :dur_infectious_non_centered, :ϕ_cases_non_centered, :ϕ_deaths_non_centered, :ρ_death_non_centered, :S_SEI_non_centered, :I_EI_non_centered),
    ESS(:R0_params_non_centered),
    ESS(:IFR_t_params_non_centered),
    ESS(:α_t_params_non_centered))

simulated_MAP_init = optimize_many_MAP(my_model_simulated, 10, 1, true)[1]

Random.seed!(seed)
MAP_noise = randn(length(MAP_init))
Random.seed!(seed)
posterior_samples = sample(my_model_simulated, alg, n_samples, discard_initial = 10_000, thin = 10, init_params = simulated_MAP_init * 0.95 + MAP_noise * 0.05)

mkpath(resultsdir("simulated_posterior_samples"))
wsave(resultsdir("simulated_posterior_samples", savename("simulated_posterior_samples", simulated_dict, "jld2")), @dict posterior_samples)
