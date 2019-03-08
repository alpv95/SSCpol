#!/bin/bash
#
#SBATCH --job-name=IC_ALP19task2
#SBATCH --output=current_%A_%a.out
#SBATCH --error=current_%A_%a.err
#SBATCH --time=540:00
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=7G
#SBATCH --partition=normal
#SBATCH --array=5-8 ##these are the SLURM task ids

module load gcc
module load python/3.6.1
module load py-numpy/1.14.3_py36
module load py-pytorch/1.0.0_py36
module load viz
module load py-matplotlib/2.1.2_py36

python3 automateSED.py $SLURM_ARRAY_TASK_ID 19 2 100 2
