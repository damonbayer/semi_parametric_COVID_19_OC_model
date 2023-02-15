# Semi-parametric modeling of SARS-CoV-2 transmission using tests, cases, deaths, and seroprevalence data

This repository contains code used in the manuscript Semi-parametric modeling of [SARS-CoV-2 transmission using tests, cases, deaths, and seroprevalence data](https://arxiv.org/abs/2009.02654) by Damon Bayer (@damonbayer), Jon Fintzi (@fintzij), Isaac Goldstein (@igoldsteinh), Volodymyr Minin (@vnminin), et al.

This model is able to account for time-varying testing policy, infection-fatality ratio, and reproduction number.
Additionally, we are able to use data from seroprevalence studies.

The main model code is available in [src/bayes_seird.jl](src/bayes_seird.jl).
