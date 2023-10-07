#!/bin/bash

#SBATCH -p free   ## run on the standard partition
#SBATCH -A bayerd ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 00:20:00   ## 1 hr run time limit
#SBATCH --mem=4G
#SBATCH -o generate_posterior_predictive_and_generated_quantities-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=201-1000

module purge
module load julia/1.8.5
cd //pub/bayerd/semi_parametric_COVID_19_OC_model/

julia --project --threads 1 scripts/simulation/oc_like/sensitivity_analysis/generate_simulated_predictive_and_generated_quantities_sensitivity.jl $SLURM_ARRAY_TASK_ID
