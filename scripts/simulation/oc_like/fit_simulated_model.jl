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

simulated_dict = @dict sim_id
simulation_results_dir(args...) = projectdir("results/simulation/oc_like", args...)
mkpath(simulation_results_dir())

simulation_data_dir(args...) = projectdir("data/simulation/oc_like", args...)
simulated_data = load(simulation_data_dir("simulated_data.jld2"))["simulated_data"][sim_id,:,1]

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

## Control Parameters
n_samples = 250
n_chains = 4

## Define Prior Constants
include(projectdir("src/prior_constants.jl"))

## Define Model
include(projectdir("src/seirdc_log_ode.jl"))

## Load Data
include(projectdir("src/load_process_data.jl"))
include(projectdir("src/bayes_seird.jl"))
data_new_deaths, data_new_cases, data_seroprev_cases = map(x -> Int.(vec(vcat(get(simulated_data, x)[x]...))), [:data_new_deaths, :data_new_cases, :data_seroprev_cases])

## Create Models
model_optimization = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-11, 1e-8)
model_sample = bayes_seird(prob, data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, Tsit5(), 1e-9, 1e-6)

MAP_init = optimize_many_MAP(model_optimization, 10, 1, true)[1]
MAP_init_values = MAP_init.values.array

alg = NUTS(-1, 0.8)

Random.seed!(sim_id)
MAP_noise = [randn(length(MAP_init_values)) for x in 1:n_chains]

Random.seed!(sim_id)
posterior_samples = sample(model_sample, alg, MCMCThreads(), n_samples, n_chains, discard_initial=5, thinning=1, init_params=repeat([MAP_init_values], n_chains) * 0.95 + MAP_noise * 0.05)

mkpath(simulation_results_dir("posterior_samples"))
wsave(simulation_results_dir("posterior_samples", savename("posterior_samples", simulated_dict, "jld2")), @dict posterior_samples)
