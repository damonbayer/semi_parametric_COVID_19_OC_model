#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 06:00:00   ## 1 hr run time limit
#SBATCH --mem=4G 
#SBATCH -o fit_simulated-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=1-200

module purge
module load julia/1.8.2
cd //dfs6/pub/bayerd/semi_parametric_COVID_19_OC_model/

julia --project --threads 1 scripts/simulation/fit_simulated_model.jl $SLURM_ARRAY_TASK_ID