# https://gist.github.com/torfjelde/7794c384d82d03c36625cd25b702b8d7

# Use packages to ensure that we trigger Requires.jl.
using Zygote: Zygote
using ReverseDiff: ReverseDiff
using ForwardDiff: ForwardDiff
using Tracker: Tracker
using Memoization: Memoization # used for ReverseDiff.jl cache.

using Turing.Core: ForwardDiffAD, ReverseDiffAD, TrackerAD, ZygoteAD, CHUNKSIZE
using BenchmarkTools
const DEFAULT_ADBACKENDS = [
    ForwardDiffAD{40}(),    # chunksize=40
    ForwardDiffAD{100}(),   # chunksize=100
    # TrackerAD(),
    # ZygoteAD(),
    ReverseDiffAD{false}(), # rdcache=false
    ReverseDiffAD{true}()   # rdcache=false
]
"""
    make_turing_suite(model; kwargs...)

Create default benchmark suite for `model`.

# Keyword arguments
- `adbackends`: a collection of adbackends to use. Defaults to `$(DEFAULT_ADBACKENDS)`.
- `run_once=true`: if `true`, the body of each benchmark will be run once to avoid
  compilation to be included in the timings (this may occur if compilation runs
  longer than the allowed time limit).
- `save_grads=false`: if `true` and `run_once` is `true`, the gradients from the initial
  execution will be saved and returned as the second return-value. This is useful if you
  want to check correctness of the gradients for different backends.
# Notes
- A separate "parameter" instance (`DynamicPPL.VarInfo`) will be created for _each test_.
  Hence if you have a particularly large model, you might want to only pass one `adbackend`
  at the time.
"""
function make_turing_suite(model; adbackends=DEFAULT_ADBACKENDS, run_once=true, save_grads=false)
    suite = BenchmarkGroup()
    suite["not_linked"] = BenchmarkGroup()
    suite["linked"] = BenchmarkGroup()

    grads = Dict(:not_linked => Dict(), :linked => Dict())

    vi_orig = DynamicPPL.VarInfo(model)
    spl = DynamicPPL.SampleFromPrior()

    for adbackend in adbackends
        println(adbackend)
        vi = DynamicPPL.VarInfo(model)
        vi[spl] = deepcopy(vi_orig[spl])

        if run_once
            ℓ, ∇ℓ = Turing.Core.gradient_logp(
                adbackend,
                vi[spl],
                vi,
                model,
                spl
            )

            if save_grads
                grads[:not_linked][adbackend] = (ℓ, ∇ℓ)
            end
        end
        suite["not_linked"]["$(adbackend)"] = @benchmarkable $(Turing.Core.gradient_logp)(
            $adbackend,
            $(vi[spl]),
            $vi,
            $model,
            $spl
        )

        # Need a separate `VarInfo` for the linked version since otherwise we risk the
        # `vi` from above being mutated.
        vi_linked = deepcopy(vi)
        DynamicPPL.link!(vi_linked, spl)
        if run_once
            ℓ, ∇ℓ = Turing.Core.gradient_logp(
                adbackend,
                vi_linked[spl],
                vi_linked,
                model,
                spl
            )
            
            if save_grads
                grads[:linked][adbackend] = (ℓ, ∇ℓ)
            end
        end
        suite["linked"]["$(adbackend)"] = @benchmarkable $(Turing.Core.gradient_logp)(
            $adbackend,
            $(vi_linked[spl]),
            $vi_linked,
            $model,
            $spl
        )

    end

    return save_grads ? (suite, grads) : suite
end

tmp = make_turing_suite(my_model)

run(tmp)
make_turing_suite(my_model_forecast)
make_turing_suite(my_model_forecast_missing)