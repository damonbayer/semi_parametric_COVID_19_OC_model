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


precision_experiment_table = CSV.read("precision_experiment_table.csv", DataFrame)

function read_group(group_id)
    model_dicts = [Dict(names(precision_experiment_table_row) .=> values(precision_experiment_table_row)) for precision_experiment_table_row in eachrow(subset(precision_experiment_table, :group_id => ByRow(x -> x == group_id)))]
    model_dicts = delete!.(model_dicts, "group_id")
    files_to_read = [resultsdir("precision_experiment", savename("posterior_samples", model_dict, "jld2")) for model_dict in model_dicts]
    posterior_samples = cat([load(file_to_read)["posterior_samples"] for file_to_read in files_to_read]..., dims=3)    
end



posterior_samples = [read_group(x) for x in 1:6]

ess = ess_rhat.(posterior_samples)

ess_sum = sum.([e.nt.ess for e in ess])
ess_per_sec_sum = sum.([e.nt.ess_per_sec for e in ess])
rhat_sum = sum.([e.nt.rhat for e in ess])

ess_sum / maximum(ess_sum)
# 1, 2, 3, 4 viable. 1 and 3 preferred
# hp model preferred

ess_per_sec_sum / maximum(ess_per_sec_sum)
# 1, 2, 3, 4 viable. 2 and 4 prefferred
# lets go hp hp
rhat_sum
# 2 and 3 preferred

sortperm(sum.([e.nt.ess_per_sec for e in ess]), rev = true)
sortperm(sum.([e.nt.ess for e in ess]), rev = true)
sortperm(sum.([e.nt.rhat for e in ess]), rev = false)


1: 4 + 2 + 5 = 11
2: 1 + 4 + 2 = 7
3: 3 + 1 + 4 = 8
4: 2 + 3 + 6 = 11
5: 6 + 6 + 3 = 15
6: 5 + 6 + 1 = 12

2
4
3
1
6
5

3
1
4
2
6
5

6
2
5
3
1
4

ess[6]
posterior_samples[6]

describe.([e.nt.rhat for e in ess])
str(ess)

describe(ess_rhat(posterior_samples).nt.rhat)

sort(sum.([e.nt.ess for e in ess]), rev = true)
max(ess[1].nt.rhat)


sortperm(median.(vec.(Array.(getindex.(posterior_samples, "lp")))), rev=true)
findmax(median.(vec.(Array.(getindex.(posterior_samples, "lp")))))
findmax(mean.(vec.(Array.(getindex.(posterior_samples, "lp")))))




sum.([e.nt.ess for e in ess])
sum.([e.nt.ess_per_sec for e in ess])
sortperm(sum.([e.nt.ess_per_sec for e in ess]), rev = true)
sortperm(sum.([e.nt.ess for e in ess]), rev = true)

