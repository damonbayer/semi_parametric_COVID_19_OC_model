#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 01:00:00   ## 1 hr run time limit
#SBATCH --mem=3G 
#SBATCH -o summarize_simulated_model_results-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=1-200

module purge
module load julia-lts
cd //dfs6/pub/bayerd/semi_parametric_COVID_19_OC_model/

julia --project --threads 1 scripts/summarize_simulated_model_results.jl $SLURM_ARRAY_TASK_ID