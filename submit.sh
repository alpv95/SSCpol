#!/bin/bash
#
#SBATCH --job-name=IC_ALP
#SBATCH --output=current_%A_%a.out
#SBATCH --error=current_%A_%a.err
#SBATCH --time=210:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=iric
#SBATCH --array=[1-3]

module load gcc
module load python/3.6.1
module load py-numpy/1.14.3_py36
module load py-pytorch/1.0.0_py36
module load viz
module load py-matplotlib/2.1.2_py36

srun python3 automateSED.py $SLURM_ARRAY_TASK_ID