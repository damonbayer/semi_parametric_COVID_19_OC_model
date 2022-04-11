chain_seed = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])

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

## Control Parameters
max_t = 42
max_t_forecast = max_t + 4
n_samples = 2000

savename_dict = (@dict chain_seed max_t max_t_forecast)

## Define Prior Constants
# This should be differentiated for different models
include("prior_constants.jl")

## Define Model
# Should be differntiated for different models
Turing.setadbackend(:forwarddiff)
include("seirdc_log_ode.jl")
include("bayes_seird.jl")

## Load Data 
include("load_process_data.jl")

my_model = bayes_seird(data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times)
my_model_forecast = bayes_seird(data_new_deaths_forecast, data_new_cases_forecast, tests_forecast, data_seroprev_cases, seroprev_tests, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast)
my_model_forecast_missing = bayes_seird(missing_new_deaths_forecast, missing_new_cases_forecast, tests_forecast, missing_seroprev_cases_forecast, seroprev_tests_forecast, obstimes_forecast, seroprev_times_forecast, param_change_times_forecast)

MAP_init = optimize_many_MAP(my_model)
Random.seed!(chain_seed)
init_noise = randn(length(MAP_init.values.array)) * 0.1

alg = Gibbs(NUTS(1000, 0.65, :dur_latent_non_centered, :dur_infectious_non_centered, :ϕ_cases_non_centered, :ϕ_deaths_non_centered, :ρ_death_non_centered, :S_SEI_non_centered, :I_EI_non_centered),
    ESS(:R0_params_non_centered),
    ESS(:IFR_t_params_non_centered),
    ESS(:α_t_params_non_centered))

Random.seed!(chain_seed)
posterior_samples = sample(my_model, alg, 50_000, discard_initial = 10000, init_params = MAP_init.values.array + init_noise)

wsave(resultsdir("posterior_samples", savename("posterior_samples", savename_dict, "jld2")), @dict posterior_samples)
# wload(resultsdir("posterior_samples", savename("posterior_samples", savename_dict, "jld2")))["posterior_samples"]

original_chains = Chains(posterior_samples, :parameters)
model = my_model
model_forecast = my_model_forecast
augment_type = "zeros"
forecast_chains_zeros = augment_chains_with_forecast_samples(Chains(posterior_samples, :parameters), my_model, my_model_forecast, "zeros")

## 
MAP_init_chains_multi =
    cat(
        map(function (i)
                Random.seed!(i)
                init_noise = randn(length(MAP_init.values.array)) * 0.1
                Chains(transpose(reshape(repeat(MAP_init.values.array + init_noise, n_samples), :, n_samples)), names(MAP_init.values, 1))
            end, 1:4)...,
        dims = 3)

MAP_init_chains = Chains(transpose(reshape(repeat(MAP_init.values.array, n_samples), :, n_samples)), names(MAP_init.values, 1))
# MAP_init_chains_multi = cat(repeat([MAP_init_chains], 4)..., dims = 3)

forecast_chains_zeros = augment_chains_with_forecast_samples(MAP_init_chains, my_model, my_model_forecast, "zeros")
forecast_chains_randn = augment_chains_with_forecast_samples(MAP_init_chains, my_model, my_model_forecast, "randn")

forecast_chains_zeros_multi = augment_chains_with_forecast_samples(MAP_init_chains_multi, my_model, my_model_forecast, "zeros")
forecast_chains_randn_multi = augment_chains_with_forecast_samples(MAP_init_chains_multi, my_model, my_model_forecast, "randn")

gq_forecast_randn = get_gq_chains(my_model_forecast, forecast_chains_randn_multi)
gq_forecast_zeros = get_gq_chains(my_model_forecast, forecast_chains_zeros_multi)


CSV.write("gq_forecast_randn.csv", DataFrame(gq_forecast_randn))
CSV.write("gq_forecast_zeros.csv", DataFrame(gq_forecast_zeros))
# Right now exploring what kind of sd we should have in the init noise
# later we will write code to save the model results
# then post processing code to read, combine, process, and re-save as csv

## Predictive Distribution
predictive_zeros = predict(my_model_forecast_missing, forecast_chains_zeros_multi)
predictive_randn = predict(my_model_forecast_missing, forecast_chains_randn_multi)

CSV.write("predictive_zeros.csv", DataFrame(predictive_zeros))
CSV.write("predictive_randn.csv", DataFrame(predictive_randn))

## Fit the model
alg = Gibbs(NUTS(1000, 0.65, :dur_latent_non_centered, :dur_infectious_non_centered, :ϕ_cases_non_centered, :ϕ_deaths_non_centered, :ρ_death_non_centered, :S_SEI_non_centered, :I_EI_non_centered),
    ESS(:R0_params_non_centered),
    ESS(:IFR_t_params_non_centered),
    ESS(:α_t_params_non_centered))

# posterior_samples = sample(my_model, alg, n_samples, thinning = 20, discard_initial = 20000, init_params = map_values + map_noise)
Random.seed!(1)
posterior_samples = sample(my_model, alg, 60_000, init_params = MAP_init.values.array + init_noise)
summarize(posterior_samples, "dur_latent_non_centered")
posterior_samples[:lp]
using StatsPlots
plot(posterior_samples[:lp])

summarize(posterior_samples[20000:100:60000, :, :], sections=[:internals])
summarize(posterior_samples[5000:1:60000, :, :], sections=[:internals])
summarize(posterior_samples[1000:1:60000, :, :], sections=[:internals])
summarize(posterior_samples[1000:25:60000, :, :], sections=[:internals])
summarize(posterior_samples, sections=[:internals])
plot(posterior_samples[20000:100:60000, :lp, :])
