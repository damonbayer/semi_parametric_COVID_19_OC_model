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

simulated_dict = @dict max_t use_tests use_seroprev constant_R0 constant_alpha constant_IFR double_IFR_0 half_alpha_0 half_S_0 half_R0_0 seed

n_samples = 10_000
n_chains = 4

## Define Prior Constants
include(projectdir("src/prior_constants.jl"))

## Define Model
Turing.setadbackend(:forwarddiff)
include(projectdir("src/seirdc_log_ode.jl"))

## Load Data
fixed_ϕ_cases_non_centered = 1.5
fixed_ϕ_deaths_non_centered = 1.5
fixed_σ_R0_non_centered = 1.5
fixed_σ_IFR_non_centered = 1.5
fixed_σ_α_non_centered = 1.5
fixed_σ_ρ_cases_non_centered = 1.5
include(projectdir("src/load_process_data.jl"))
include(projectdir("src/bayes_seird_fixed_phi.jl"))

## Create Models
# I don't remember why we needed phixed phi models in the first place.
# Maybe because the MAP estimates were weird?
my_model_fixed_phi = bayes_seird_fixed_phi(data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, false)
MAP_init_fixed_phi = optimize_many_MAP(my_model_fixed_phi, 10, 1, true)[1]
named_MAP_init_fixed_phi = (;zip(flattened_varnames_list(my_model_fixed_phi), MAP_init_fixed_phi)...)

named_MAP_init_unordered = merge(named_MAP_init_fixed_phi,
    (;zip((:ϕ_cases_non_centered, :ϕ_deaths_non_centered, Symbol("R0_params_non_centered[2]"), Symbol("IFR_t_params_non_centered[2]"), Symbol("α_t_params_non_centered[2]"), Symbol("ρ_cases_t_params_non_centered[2]")),
    [fixed_ϕ_cases_non_centered, fixed_ϕ_deaths_non_centered, fixed_σ_R0_non_centered, fixed_σ_IFR_non_centered, fixed_σ_α_non_centered, fixed_σ_ρ_cases_non_centered])...))

include(projectdir("src/bayes_seird.jl"))

my_model = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, false)
my_model_gq = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, true)
my_model_forecast_missing = bayes_seird(prob, missing_new_deaths_forecast, missing_new_cases_forecast, tests_forecast, missing_seroprev_cases_forecast, seroprev_tests_forecast, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, true)

true_parameters = [named_MAP_init_unordered[param] for param in flattened_varnames_list(my_model)]
# Save true parameters as jld2
wsave(datadir("simulated_data", savename("true_parameters", simulated_dict, "jld2")), @dict true_parameters)

true_parameters_chain = Chains([true_parameters], flattened_varnames_list(my_model))
# Save true parameters as csv
CSV.write(datadir("simulated_data", savename("true_parameters", simulated_dict, "csv")), DataFrame(true_parameters_chain))

true_generated_quantities_chain = get_gq_chains(my_model_gq, true_parameters_chain)
CSV.write(datadir("simulated_data", savename("true_generated_quantities", simulated_dict, "csv")), DataFrame(true_generated_quantities_chain))

true_parameters_chain_for_data_gen = Chains(transpose(hcat(repeat([true_parameters], 200)...)), flattened_varnames_list(my_model))
Random.seed!(1)
data_chain = predict(my_model_forecast_missing, true_parameters_chain_for_data_gen)
# Save simulated data
wsave(datadir("simulated_data", savename("simulated_data", simulated_dict, "jld2")), @dict data_chain)
CSV.write(datadir("simulated_data", savename("simulated_data", simulated_dict, "csv")), DataFrame(data_chain))

Random.seed!(1)
prior_samples = sample(my_model, Prior(), MCMCThreads(), n_samples, n_chains)
prior_samples_summary = innerjoin(DataFrame.(describe(prior_samples, q = [0.1, 0.9]))..., on = :parameters)
# Save prior samples summary
CSV.write(resultsdir("simulated_posterior_samples_summary", savename("simulated_prior_samples_summary", simulated_dict, "csv")), prior_samples_summary)


generated_quantities_prior = get_gq_chains(my_model_gq, prior_samples)
generated_quantities_prior_summary = innerjoin(DataFrame.(describe(generated_quantities_prior, q = [0.1, 0.9]))..., on = :parameters)
# Save prior generated_quantities summary
CSV.write(resultsdir("simulated_generated_quantities_summary", savename("simulated_prior_generated_quantities_summary", simulated_dict, "csv")), generated_quantities_prior_summary)
