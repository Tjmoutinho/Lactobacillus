#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 10:00:00
#SBATCH --output=error_log.out
#SBATCH --mem=50000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard

module load gcc openmpi anaconda2/5.1.0 gurobi/8.0.0

python test_re_gapfill.py
