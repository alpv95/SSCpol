#!/bin/bash
#
#SBATCH --job-name=7final6
#SBATCH --output=current_%A_%a.out
#SBATCH --error=current_%A_%a.err
#SBATCH --time=400:00
##SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9G
#SBATCH --partition=iric
#SBATCH --array=1-200 ##these are the SLURM task ids

module load gcc
module load python/3.6.1
module load py-numpy/1.14.3_py36
module load py-pytorch/1.0.0_py36
module load viz
module load py-matplotlib/2.1.2_py36

#gcc jet_synconlyALPWFlux.c mtwister.c mtwister.h jet_fns.c jet_fns.h -lm
python3 automateSED.py $SLURM_ARRAY_TASK_ID 7 1 90 1 6.0

