model_id = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])

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

model_dict = subset(CSV.read("model_table.csv", DataFrame), :model_id => ByRow(x -> x == model_id))[1,:]|> x -> Dict(names(x) .=> values(x))

max_t = float(model_dict["max_t"])
seed = model_dict["seed"]
use_tests = model_dict["use_tests"]
use_seroprev = model_dict["use_seroprev"]
constant_R0 = model_dict["constant_R0"]
constant_alpha = model_dict["constant_alpha"]
constant_IFR = model_dict["constant_IFR"]
double_IFR_0 = model_dict["double_IFR_0"]
half_alpha_0 = model_dict["half_alpha_0"]
half_S_0 = model_dict["half_S_0"]
half_R0_0 = model_dict["half_R0_0"]

## Control Parameters
# max_t_forecast = max_t + 4
n_samples = 10_000

## Define Prior Constants
# This should be differentiated for different models
include(projectdir("src/prior_constants.jl"))

## Define Model
# Should be differntiated for different models
Turing.setadbackend(:forwarddiff)
include(projectdir("src/seirdc_log_ode.jl"))

## Load Data 
include(projectdir("src/load_process_data.jl"))
include(projectdir("src/bayes_seird.jl"))

## Create Models
my_model = bayes_seird(data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, false)
# my_model_forecast = bayes_seird(data_new_deaths_forecast, data_new_cases_forecast, tests_forecast, data_seroprev_cases, seroprev_tests, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, true)
# my_model_forecast_missing = bayes_seird(missing_new_deaths_forecast, missing_new_cases_forecast, tests_forecast, missing_seroprev_cases_forecast, seroprev_tests_forecast, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast, use_tests, true, constant_R0, constant_alpha, constant_IFR, true)
MAP_init = optimize_many_MAP(my_model, 100, 4, true)

if use_tests
    alg = Gibbs(NUTS(1000, 0.65, :dur_latent_non_centered, :dur_infectious_non_centered, :ϕ_cases_non_centered, :ϕ_deaths_non_centered, :ρ_death_non_centered, :S_SEI_non_centered, :I_EI_non_centered),
        ESS(:R0_params_non_centered),
        ESS(:IFR_t_params_non_centered),
        ESS(:α_t_params_non_centered))
else
    alg = Gibbs(NUTS(1000, 0.65, :dur_latent_non_centered, :dur_infectious_non_centered, :ϕ_cases_non_centered, :ϕ_deaths_non_centered, :ρ_death_non_centered, :S_SEI_non_centered, :I_EI_non_centered),
        ESS(:R0_params_non_centered),
        ESS(:IFR_t_params_non_centered),
        ESS(:ρ_cases_t_params_non_centered))
end

Random.seed!(seed)
posterior_samples = sample(my_model, alg, n_samples, discard_initial = 10_000, thin = 5, init_params = MAP_init[seed] * 0.95)

wsave(resultsdir("posterior_samples", savename("posterior_samples", model_dict, "jld2")), @dict posterior_samples)
