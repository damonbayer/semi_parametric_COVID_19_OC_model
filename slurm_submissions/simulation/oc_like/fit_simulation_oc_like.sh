#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 16:00:00   ## 1 hr run time limit
#SBATCH --mem=16G
#SBATCH -o fit_simulation_oc_like-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=201-1000

module purge
module load julia/1.8.5
cd //pub/bayerd/semi_parametric_COVID_19_OC_model/

julia --project --threads 4 scripts/simulation/oc_like/sensitivity_analysis/fit_simulated_model_sensitivity.jl $SLURM_ARRAY_TASK_ID
