# Semi-parametric modeling of SARS-CoV-2 transmission using tests, cases, deaths, and seroprevalence data

This repository contains code used in the manuscript Semi-parametric modeling of [SARS-CoV-2 transmission using tests, cases, deaths, and seroprevalence data](https://arxiv.org/abs/2009.02654) by Damon Bayer ([@damonbayer](https://github.com/damonbayer)), Jon Fintzi ([@fintzij](https://github.com/fintzij)), Isaac Goldstein ([@igoldsteinh](https://github.com/igoldsteinh)), Volodymyr Minin ([@vnminin](https://github.com/vnminin)), et al.

This model is able to account for time-varying testing policy, infection-fatality ratio, and reproduction number.
Additionally, we are able to use data from seroprevalence studies.

The main model code is available in [src/bayes_seird.jl](src/bayes_seird.jl).

The workflow for fitting and evaluating models is available in the [slurm_submissions](slurm_submissions) directory.

In order, we run
- [fit_models.sh](slurm_submissions/fit_models.sh)
- [generate_predictive_and_generated_quantities.sh](slurm_submissions/generate_predictive_and_generated_quantities.sh)
- [tidy_predictive_and_generated_quantities.sh](slurm_submissions/tidy_predictive_and_generated_quantities.sh)
