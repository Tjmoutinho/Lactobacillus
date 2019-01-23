#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 100:00:00
#SBATCH --output=error_log.out
#SBATCH --mem=12000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard

module load bioconda/py2.7 gcc openmpi gurobi/8.0.1

python Create_likelihood_files.py
