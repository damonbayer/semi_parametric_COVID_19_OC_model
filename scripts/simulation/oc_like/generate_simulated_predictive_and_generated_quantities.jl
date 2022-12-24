sim_id = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])

using Revise
using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using LogExpFunctions
using Turing
using LinearAlgebra
using FillArrays
using Random
using JLD2
using FileIO
using DrWatson
using semi_parametric_COVID_19_OC_model

simulation_data_dir(args...) = projectdir("data/simulation/oc_like", args...)
simulation_results_dir(args...) = projectdir("results/simulation/oc_like", args...)

simulated_data = load(simulation_data_dir("simulated_data.jld2"))["simulated_data"][sim_id,:,1]

mkpath(simulation_results_dir("posterior_predictive"))
mkpath(simulation_results_dir("posterior_generated_quantities"))
mkpath(simulation_results_dir("duration"))

simulated_dict = @dict sim_id

posterior_samples = load(simulation_results_dir("posterior_samples", savename("posterior_samples", simulated_dict, "jld2")))["posterior_samples"]

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

durations = DataFrame(wall=MCMCChains.wall_duration(posterior_samples), compute = MCMCChains.compute_duration(posterior_samples))
CSV.write(simulation_results_dir("duration", savename("duration", simulated_dict, "csv")), durations)

## Define Prior Constants
include(projectdir("src/prior_constants.jl"))

## Define Model
include(projectdir("src/seirdc_log_ode.jl"))

## Load Data
include(projectdir("src/load_process_data.jl"))
include(projectdir("src/bayes_seird.jl"))
data_new_deaths, data_new_cases, data_seroprev_cases = map(x -> Int.(vec(vcat(get(simulated_data, x)[x]...))), [:data_new_deaths, :data_new_cases, :data_seroprev_cases])

model_sample = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-9, 1e-6)
model_generated_quantities = bayes_seird(prob, data_new_deaths_forecast, data_new_cases_forecast, tests_forecast, data_seroprev_cases, seroprev_tests, obstimes_forecast, seroprev_times, param_change_times_forecast, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-11, 1e-8)
model_predict = bayes_seird(prob, missing_new_deaths_forecast, missing_new_cases_forecast, tests_forecast, missing_seroprev_cases_forecast, seroprev_tests, obstimes_forecast, seroprev_times, param_change_times_forecast, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-11, 1e-8)

## Augment Samples for forecasting
augmented_posterior_samples = augment_chains_with_forecast_samples(posterior_samples, model_sample, model_generated_quantities, "randn")

Random.seed!(1)
posterior_predictive = predict(model_predict, augmented_posterior_samples)
CSV.write(simulation_results_dir("posterior_predictive", savename("posterior_predictive", simulated_dict, "csv")), DataFrame(posterior_predictive))

posterior_generated_quantities = hcat(Chains(generated_quantities(model_generated_quantities, augmented_posterior_samples)), getlogp_chain(posterior_samples))
CSV.write(simulation_results_dir("posterior_generated_quantities", savename("posterior_generated_quantities", simulated_dict, "csv")), DataFrame(posterior_generated_quantities))
