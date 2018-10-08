#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 10:00:00
#SBATCH --array=0-1165
#SBATCH --output=error_logs/error_log_%a.out
#SBATCH --mem=6000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard

module load gcc openmpi anaconda2/5.1.0

python pFBA_likelihoods_array.py /scratch/tjm4k/Lactobacillus/array_input_V3/input_${SLURM_ARRAY_TASK_ID}.txt
