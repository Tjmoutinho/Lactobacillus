#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 100:00:00
#SBATCH --array=0-179
#SBATCH --output=error_logs/error_log_%a.out
#SBATCH --mem=12000
#SBATCH --account=CSBLRivanna
#SBATCH --partition=standard

module load gcc openmpi anaconda2/5.1.0 

python gapfill_array.py /scratch/tjm4k/Lactobacillus/array_input/input_${SLURM_ARRAY_TASK_ID}.txt
