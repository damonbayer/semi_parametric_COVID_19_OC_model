# Simulate data to look like OC data
seed = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])
model_design = 46

using Revise
using CSV
using DataFrames
using DifferentialEquations
using LogExpFunctions
using Turing
using LinearAlgebra
using FillArrays
using Random
using JLD2
using FileIO
using DrWatson
using semi_parametric_COVID_19_OC_model

simulation_results_dir(args...) = projectdir("results/simulation/oc_like", args...)
mkpath(simulation_results_dir())
simulation_data_dir(args...) = projectdir("data/simulation/oc_like", args...)
mkpath(simulation_data_dir())

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

# We will simulate priors once since they are the same for every simulation
n_samples = 2000
n_datasets = 200

## Define Prior Constants
include(projectdir("src/prior_constants.jl"))

## Define Model
include(projectdir("src/seirdc_log_ode.jl"))

## Load Data
include(projectdir("src/load_process_data.jl"))
include(projectdir("src/bayes_seird.jl"))

model_optimization = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-11, 1e-8)
model_sample = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-9, 1e-6)
model_generated_quantities = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-11, 1e-8)
model_predict = bayes_seird(prob, missing_new_deaths_forecast, missing_new_cases_forecast, tests, missing_seroprev_cases_forecast, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-11, 1e-8)

## Save Priors
Random.seed!(seed)
prior_samples = sample(model_sample, Prior(), n_samples)
mkpath(simulation_results_dir("prior_samples"))
CSV.write(simulation_results_dir("prior_samples/prior_samples.csv"), DataFrame(prior_samples))

Random.seed!(seed)
prior_generated_quantities = Chains(generated_quantities(model_generated_quantities, prior_samples))
mkpath(simulation_results_dir("prior_generated_quantities"))
CSV.write(simulation_results_dir("prior_generated_quantities/prior_generated_quantities.csv"), DataFrame(prior_generated_quantities))

Random.seed!(seed)
prior_predictive = predict(model_predict, prior_samples)
mkpath(simulation_results_dir("prior_predictive"))
CSV.write(simulation_results_dir("prior_predictive/prior_predictive.csv"), DataFrame(prior_predictive))

## Optimize MAP
# MAP_init_sim = optimize_many_MAP(model_optimization, 10, 1, true)
# MAP_init_names_sim = names(MAP_init_sim[1].values)[1]
# MAP_init_values_sim = MAP_init_sim[1].values.array

# ## Save true parameters
# true_parameters_chain = Chains(transpose(hcat(repeat([MAP_init_values_sim], 1)...)), MAP_init_names_sim);
# CSV.write(simulation_data_dir("true_parameters_chain.csv"), DataFrame(true_parameters_chain))

# true_generated_quantities_chain = Chains(generated_quantities(model_generated_quantities, true_parameters_chain))
# CSV.write(simulation_data_dir("true_generated_quantities.csv"), DataFrame(true_generated_quantities_chain))

# ## Simulate and save data
# true_parameters_chains_for_predict = Chains(transpose(hcat(repeat([MAP_init_values_sim], n_datasets)...)), MAP_init_names_sim);
# Random.seed!(seed)
# simulated_data = predict(model_predict, true_parameters_chains_for_predict)
# CSV.write(simulation_data_dir("simulated_data.csv"), DataFrame(simulated_data))
# wsave(simulation_data_dir("simulated_data.jld2"), @dict simulated_data)

model_table = CSV.read("model_table.csv", DataFrame)
model_dicts = [Dict(names(model_table_row) .=> values(model_table_row)) for model_table_row in eachrow(subset(model_table, :model_design => ByRow(x -> x == model_design)))]

posterior_files_to_read = [resultsdir("posterior_samples", savename("posterior_samples", model_dict, "jld2")) for model_dict in model_dicts]
posterior_samples = cat([load(file_to_read)["posterior_samples"] for file_to_read in posterior_files_to_read]..., dims=3)
posterior_summary = summarize(posterior_samples)

## Save true parameters
true_parameters_chain = Chains(transpose(hcat(repeat([posterior_summary.nt.mean], 1)...)), posterior_summary.nt.parameters);
CSV.write(simulation_data_dir("true_parameters_chain.csv"), DataFrame(true_parameters_chain))

true_generated_quantities_chain = Chains(generated_quantities(model_generated_quantities, true_parameters_chain))
CSV.write(simulation_data_dir("true_generated_quantities.csv"), DataFrame(true_generated_quantities_chain))

## Simulate and save data
true_parameters_chains_for_predict = Chains(transpose(hcat(repeat([posterior_summary.nt.mean], n_datasets)...)), posterior_summary.nt.parameters);
Random.seed!(seed)
simulated_data = predict(model_predict, true_parameters_chains_for_predict)
CSV.write(simulation_data_dir("simulated_data.csv"), DataFrame(simulated_data))
wsave(simulation_data_dir("simulated_data.jld2"), @dict simulated_data)