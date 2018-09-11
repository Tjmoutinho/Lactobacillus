#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 10:00:00
#SBATCH --output=error_log.out
#SBATCH --mem=120000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=gpu
#SBATCH --gres=gpu:p100:1

module load gcc openmpi anaconda2/5.1.0 gurobi/8.0.0

python probanno_gapfill_10_models.py
