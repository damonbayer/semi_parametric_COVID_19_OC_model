seed = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])

using Revise
using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using LogExpFunctions
using Turing
using LinearAlgebra
using Random
using JLD2
using FileIO
using DrWatson
using semi_parametric_COVID_19_OC_model

mkpath(resultsdir("simulated_posterior_samples"))

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
seed = 1
simulated_dict_seed_1 = @dict seed max_t use_tests use_seroprev constant_R0 constant_alpha constant_IFR double_IFR_0 half_alpha_0 half_S_0 half_R0_0
seed = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])

## Control Parameters
n_samples = 10_000
n_chains = 4

## Define Prior Constants
include(projectdir("src/prior_constants.jl"))

## Define Model
Turing.setadbackend(:forwarddiff)
include(projectdir("src/seirdc_log_ode.jl"))

## Load Data
include(projectdir("src/load_process_data.jl"))
include(projectdir("src/bayes_seird.jl"))

## Create Models
my_model = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, false)
my_model_forecast_missing = bayes_seird(prob, missing_new_deaths_forecast, missing_new_cases_forecast, tests_forecast, missing_seroprev_cases_forecast, seroprev_tests_forecast, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, true)
data_chain = load(datadir("simulated_data", savename("simulated_data", simulated_dict_seed_1, "jld2")))["data_chain"][seed]

simulated_data_new_deaths = round.(Int, vcat(get(data_chain, :data_new_deaths)[:data_new_deaths]...))
simulated_data_new_cases = round.(Int, vcat(get(data_chain, :data_new_cases)[:data_new_cases]...))
simulated_data_seroprev_cases = round.(Int, vcat(get(data_chain, :data_seroprev_cases)[:data_seroprev_cases]...))

my_model_simulated = bayes_seird(prob, simulated_data_new_deaths, simulated_data_new_cases, tests_forecast, simulated_data_seroprev_cases, seroprev_tests_forecast, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, false)
my_model_simulated_gq = bayes_seird(prob, simulated_data_new_deaths, simulated_data_new_cases, tests_forecast, simulated_data_seroprev_cases, seroprev_tests_forecast, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, true)

alg = Gibbs(NUTS(-1, 0.65, :dur_latent_non_centered, :dur_infectious_non_centered, :ϕ_cases_non_centered, :ϕ_deaths_non_centered, :ρ_death_non_centered, :S_SEI_non_centered, :I_EI_non_centered),
    ESS(:R0_params_non_centered),
    ESS(:IFR_t_params_non_centered),
    ESS(:α_t_params_non_centered))

MAP_init = optimize_many_MAP(my_model_simulated, 10, 1, true)[1]

Random.seed!(seed)
MAP_noise = [randn(length(MAP_init)) for x in 1:n_chains]

Random.seed!(seed)
posterior_samples = sample(my_model_simulated, alg, MCMCThreads(), n_samples, n_chains, discard_initial = 10_000, thin = 10, init_params = repeat([MAP_init], n_chains) * 0.95 + MAP_noise * 0.05)
posterior_samples_summary = innerjoin(DataFrame.(describe(posterior_samples, q = [0.1, 0.5, 0.9]))..., on = :parameters)
CSV.write(resultsdir("simulated_posterior_samples_summary", savename("simulated_posterior_samples_summary", simulated_dict, "csv")), posterior_samples_summary)

generated_quantities = get_gq_chains(my_model_simulated_gq, posterior_samples)
generated_quantities_summary = innerjoin(DataFrame.(describe(generated_quantities, q = [0.1, 0.5, 0.9]))..., on = :parameters)
CSV.write(resultsdir("simulated_generated_quantities_summary", savename("simulated_generated_quantities_summary", simulated_dict, "csv")), generated_quantities_summary)