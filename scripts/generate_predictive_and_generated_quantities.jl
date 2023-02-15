# Generate predictive and generated quantities for a model with a given model_design
model_design = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])

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

mkpath(resultsdir("prior_predictive"))
mkpath(resultsdir("prior_generated_quantities"))
mkpath(resultsdir("posterior_predictive"))
mkpath(resultsdir("posterior_generated_quantities"))
mkpath(resultsdir("duration"))

model_table = CSV.read("model_table.csv", DataFrame)

model_dicts = [Dict(names(model_table_row) .=> values(model_table_row)) for model_table_row in eachrow(subset(model_table, :model_design => ByRow(x -> x == model_design)))]

posterior_files_to_read = [resultsdir("posterior_samples", savename("posterior_samples", model_dict, "jld2")) for model_dict in model_dicts]
posterior_samples = cat([load(file_to_read)["posterior_samples"] for file_to_read in posterior_files_to_read]..., dims=3)

prior_files_to_read = [resultsdir("prior_samples", savename("prior_samples", model_dict, "jld2")) for model_dict in model_dicts]
prior_samples = cat([load(file_to_read)["prior_samples"] for file_to_read in prior_files_to_read]..., dims=3)

model_dict = model_dicts[1]
delete!(model_dict, "seed")
max_t = float(model_dict["max_t"])
use_tests = model_dict["use_tests"]
use_seroprev = model_dict["use_seroprev"]
constant_R0 = model_dict["constant_R0"]
constant_alpha = model_dict["constant_alpha"]
constant_IFR = model_dict["constant_IFR"]
double_IFR_0 = model_dict["double_IFR_0"]
half_alpha_0 = model_dict["half_alpha_0"]
half_S_0 = model_dict["half_S_0"]
half_R0_0 = model_dict["half_R0_0"]

durations = DataFrame(wall=MCMCChains.wall_duration(posterior_samples), compute = MCMCChains.compute_duration(posterior_samples))
CSV.write(resultsdir("duration", savename("duration", model_dict, "csv")), durations)

max_t_forecast = max_t + 4

include(projectdir("src/prior_constants.jl"))
include(projectdir("src/seirdc_log_ode.jl"))
include(projectdir("src/load_process_data.jl"))
include(projectdir("src/bayes_seird.jl"))

model_sample = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-9, 1e-6)
model_generated_quantities = bayes_seird(prob, data_new_deaths_forecast, data_new_cases_forecast, tests_forecast, data_seroprev_cases, seroprev_tests, obstimes_forecast, seroprev_times, param_change_times_forecast, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-11, 1e-8)
model_predict = bayes_seird(prob, missing_new_deaths_forecast, missing_new_cases_forecast, tests_forecast, missing_seroprev_cases_forecast, seroprev_tests, obstimes_forecast, seroprev_times, param_change_times_forecast, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-11, 1e-8)

## Augment Samples for forecasting
augmented_prior_samples = augment_chains_with_forecast_samples(prior_samples, model_sample, model_generated_quantities, "randn")
augmented_posterior_samples = augment_chains_with_forecast_samples(posterior_samples, model_sample, model_generated_quantities, "randn")

Random.seed!(1)
prior_predictive = predict(model_predict, augmented_prior_samples)
CSV.write(resultsdir("prior_predictive", savename("prior_predictive", model_dict, "csv")), DataFrame(prior_predictive))

prior_generated_quantities = hcat(Chains(generated_quantities(model_generated_quantities, augmented_prior_samples)), getlogp_chain(prior_samples))
CSV.write(resultsdir("prior_generated_quantities", savename("prior_generated_quantities", model_dict, "csv")), DataFrame(prior_generated_quantities))

Random.seed!(1)
posterior_predictive = predict(model_predict, augmented_posterior_samples)
CSV.write(resultsdir("posterior_predictive", savename("posterior_predictive", model_dict, "csv")), DataFrame(posterior_predictive))

posterior_generated_quantities = hcat(Chains(generated_quantities(model_generated_quantities, augmented_posterior_samples)), getlogp_chain(posterior_samples))
CSV.write(resultsdir("posterior_generated_quantities", savename("posterior_generated_quantities", model_dict, "csv")), DataFrame(posterior_generated_quantities))
