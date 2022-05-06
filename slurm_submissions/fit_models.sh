#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 06:00:00   ## 1 hr run time limit
#SBATCH --mem=4G 
#SBATCH -o fit_models-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375

module purge
module load julia-lts
cd //dfs6/pub/bayerd/semi_parametric_COVID_19_OC_model/

julia --project --threads 1 scripts/fit_model.jl $SLURM_ARRAY_TASK_ID