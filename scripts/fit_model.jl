# Fit a model with a given model_design
model_id = isempty(ARGS) ? 309 : parse(Int64, ARGS[1])

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

model_dict = subset(CSV.read("model_table.csv", DataFrame), :model_id => ByRow(x -> x == model_id))[1, :] |> x -> Dict(names(x) .=> values(x))

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
n_samples = 250

## Define Prior Constants
include(projectdir("src/prior_constants.jl"))

## Define Model
include(projectdir("src/seirdc_log_ode.jl"))

## Load Data
include(projectdir("src/load_process_data.jl"))
include(projectdir("src/bayes_seird.jl"))

## Create Models
model_optimization = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-11, 1e-8)
model_sample = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-9, 1e-6)

# Sample Prior
Random.seed!(seed)
prior_samples = sample(model_sample, Prior(), n_samples)
mkpath(resultsdir("prior_samples"))
wsave(resultsdir("prior_samples", savename("prior_samples", model_dict, "jld2")), @dict prior_samples)

# Fit Posterior
MAP_init = optimize_many_MAP(model_optimization, 10, 1, true)[1]
MAP_init_values = MAP_init.values.array

alg = NUTS(-1, 0.8)

Random.seed!(seed)
MAP_noise = randn(length(MAP_init_values))
Random.seed!(seed)
posterior_samples = sample(model_sample, alg, n_samples, discard_initial=5, thinning=1, init_params=MAP_init_values * 0.95 + MAP_noise * 0.05)
mkpath(resultsdir("posterior_samples"))
wsave(resultsdir("posterior_samples", savename("posterior_samples", model_dict, "jld2")), @dict posterior_samples)
