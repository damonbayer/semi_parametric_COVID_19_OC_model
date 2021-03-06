#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 1          ## request 4 tasks (4 CPUs)
#SBATCH -t 04:00:00   ## 1 hr run time limit
#SBATCH --mem=11G
#SBATCH -o tidy_posterior_predictive-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu
#SBATCH --array=113-135

module purge
module load R
cd //dfs6/pub/bayerd/semi_parametric_COVID_19_OC_model/

Rscript scripts/tidy_posterior_predictive.R $SLURM_ARRAY_TASK_ID
