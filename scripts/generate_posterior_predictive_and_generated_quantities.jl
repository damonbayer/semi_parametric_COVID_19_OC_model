model_design = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])
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

mkpath(resultsdir("posterior_predictive"))
mkpath(resultsdir("generated_quantities"))

model_table = CSV.read("model_table.csv", DataFrame)

model_dicts = [Dict(names(model_table_row) .=> values(model_table_row)) for model_table_row in eachrow(subset(model_table, :model_design => ByRow(x -> x == model_design)))]
files_to_read = [resultsdir("posterior_samples", savename("posterior_samples", model_dict, "jld2")) for model_dict in model_dicts]
posterior_samples = cat([load(file_to_read)["posterior_samples"] for file_to_read in files_to_read]..., dims=3)

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

max_t_forecast = max_t + 4

include(projectdir("src/prior_constants.jl"))
include(projectdir("src/seirdc_log_ode.jl"))
include(projectdir("src/load_process_data.jl"))
include(projectdir("src/bayes_seird.jl"))

# Let actually save and load these from fit_model so we know the data is the same?
my_model = bayes_seird(data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, false)
my_model_forecast = bayes_seird(data_new_deaths_forecast, data_new_cases_forecast, tests_forecast, data_seroprev_cases, seroprev_tests, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, true)
my_model_forecast_missing = bayes_seird(missing_new_deaths_forecast, missing_new_cases_forecast, tests_forecast, missing_seroprev_cases_forecast, seroprev_tests_forecast, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, true)

augmented_posterior_samples = augment_chains_with_forecast_samples(posterior_samples, my_model, my_model_forecast, "randn")

Random.seed!(1)
posterior_predictive = predict(my_model_forecast_missing, augmented_posterior_samples)
CSV.write(resultsdir("posterior_predictive", savename("posterior_predictive", model_dict, "csv")), DataFrame(posterior_predictive))

Random.seed!(1)
generated_quantities = get_gq_chains(my_model_forecast, augmented_posterior_samples)
CSV.write(resultsdir("generated_quantities", savename("generated_quantities", model_dict, "csv")), DataFrame(generated_quantities))