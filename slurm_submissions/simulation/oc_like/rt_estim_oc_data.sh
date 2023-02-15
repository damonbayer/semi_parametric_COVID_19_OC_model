#!/bin/bash

#SBATCH -p standard   ## run on the standard partition
#SBATCH -A vminin_lab ## account to charge
#SBATCH -N 1          ## run on a single node
#SBATCH -n 4          ## request 4 tasks (4 CPUs)
#SBATCH -t 01:00:00   ## 1 hr run time limit
#SBATCH --mem=4G
#SBATCH -o rt_estim_oc_data-%A-%a.out
#SBATCH --mail-type=begin,end
#SBATCH --mail-user=bayerd@uci.edu

module purge
module load R
cd //dfs6/pub/bayerd/semi_parametric_COVID_19_OC_model/

Rscript scripts/rt_comparison_oc_data.R
