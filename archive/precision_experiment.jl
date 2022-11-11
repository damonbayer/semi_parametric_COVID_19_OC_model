experiment_id = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])

using Revise
using CSV
using DataFrames
using DifferentialEquations
using DiffEqCallbacks
using LogExpFunctions
using LinearAlgebra
using Turing
using Random
using JLD2
using FileIO
using DrWatson
using semi_parametric_COVID_19_OC_model

experiment_dict = subset(CSV.read("precision_experiment_table.csv", DataFrame), :experiment_id => ByRow(x -> x == experiment_id))[1,:]|> x -> Dict(names(x) .=> values(x))

init_type = experiment_dict["init_type"]
model_type = experiment_dict["model_type"]
seed = experiment_dict["seed"]
alg_type = experiment_dict["alg_type"]

max_t = 42.0
use_tests = true
use_seroprev = true
constant_R0 = false
constant_alpha = false
constant_IFR = false
double_IFR_0 = false
half_alpha_0 = false
half_S_0 = false
half_R0_0 =false

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
if model_type == "lp"
    my_model = bayes_seird(data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, false)
elseif model_type == "hp"
    my_model = bayes_seird(data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, true)
end

if init_type == "lp"
    MAP_init = optimize_many_MAP(my_model, 100, 1, true)
elseif init_type == "hp"
    MAP_init = optimize_many_MAP(my_model, 10, 1, true)
end

if alg_type == "Gibbs"
    alg = Gibbs(NUTS(1000, 0.65, :dur_latent_non_centered, :dur_infectious_non_centered, :ϕ_cases_non_centered, :ϕ_deaths_non_centered, :ρ_death_non_centered, :S_SEI_non_centered, :I_EI_non_centered),
    ESS(:R0_params_non_centered),
    ESS(:IFR_t_params_non_centered),
    ESS(:α_t_params_non_centered))
elseif  alg_type == "NUTS"
    alg = NUTS()
end

if init_type ∈ ["lp", "hp"]
    Random.seed!(seed)
    MAP_noise = randn(length(MAP_init[1]))
    Random.seed!(seed)
    posterior_samples = sample(my_model, alg, n_samples, discard_initial = 10_000, thinning = 1, init_params = MAP_init[1] * 0.95 + MAP_noise * 0.05)
elseif init_type == "rand"
    Random.seed!(seed)
    posterior_samples = sample(my_model, alg, n_samples, discard_initial = 10_000, thinning = 1)
end

mkpath(resultsdir("precision_experiement"))
wsave(resultsdir("precision_experiement", savename("posterior_samples", experiment_dict, "jld2")), @dict posterior_samples)
