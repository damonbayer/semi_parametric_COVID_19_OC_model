seed = isempty(ARGS) ? 1 : parse(Int64, ARGS[1])

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

mkpath(resultsdir("simulated_posterior_samples_summary"))
mkpath(resultsdir("simulated_generated_quantities_summary"))

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

simulated_dict = @dict seed max_t use_tests use_seroprev constant_R0 constant_alpha constant_IFR double_IFR_0 half_alpha_0 half_S_0 half_R0_0

if seed == 1
    prior_samples_path = resultsdir("simulated_posterior_samples", savename("simulated_prior_samples", simulated_dict, "jld2"))
    prior_samples = load(prior_samples_path)["prior_samples"]
    prior_samples_summary = innerjoin(DataFrame.(describe(prior_samples, q = [0.1, 0.9]))..., on = :parameters)
    CSV.write(replace(prior_samples_path, "simulated_prior_samples"=>"simulated_prior_samples_summary") |> x -> replace(x, "jld2"=>"csv"), prior_samples_summary)
end

posterior_samples_path = resultsdir("simulated_posterior_samples", savename("simulated_posterior_samples", simulated_dict, "jld2"))
posterior_samples = load(posterior_samples_path)["posterior_samples"]

posterior_samples_summary = innerjoin(DataFrame.(describe(posterior_samples, q = [0.1, 0.9]))..., on = :parameters)
CSV.write(replace(posterior_samples_path, "simulated_posterior_samples"=>"simulated_posterior_samples_summary") |> x -> replace(x, "jld2"=>"csv"), posterior_samples_summary)

data = CSV.read(projectdir("data", "simulated_data", savename("simulated_data", simulated_dict, "csv")), DataFrame)

include(projectdir("src/prior_constants.jl"))
include(projectdir("src/seirdc_log_ode.jl"))
include(projectdir("src/load_process_data.jl"))
data_new_deaths = round.(Int, [data[1, "data_new_deaths[" * string(i) * "]"] for i in 1:round(Int, max_t)])
data_new_cases = round.(Int, [data[1, "data_new_cases[" * string(i) * "]"] for i in 1:round(Int, max_t)])
data_seroprev_cases = round.(Int, [data[1, "data_seroprev_cases[" * string(i) * "]"] for i in 1])
include(projectdir("src/bayes_seird.jl"))

my_model = bayes_seird(data_new_deaths, data_new_cases, tests, data_seroprev_cases, seroprev_tests, obstimes, seroprev_times, param_change_times, use_tests, use_seroprev, constant_R0, constant_alpha, constant_IFR, true)

if seed == 1
    Random.seed!(1)
    generated_quantities_prior = get_gq_chains(my_model, prior_samples)
    CSV.write(replace(prior_samples_path, "simulated_prior_samples"=>"simulated_prior_generated_quantities_summary") |> x -> replace(x, "jld2"=>"csv"), generated_quantities_prior)
end

Random.seed!(1)
generated_quantities = get_gq_chains(my_model, posterior_samples)

generated_quantities_summary = innerjoin(DataFrame.(describe(generated_quantities, q = [0.1, 0.9]))..., on = :parameters)
CSV.write(replace(posterior_samples_path, "simulated_posterior_samples"=>"simulated_generated_quantities_summary") |> x -> replace(x, "jld2"=>"csv"), generated_quantities_summary)