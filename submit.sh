#!/bin/bash
#
#SBATCH --job-name=S5
#SBATCH --output=slurm/%A.out
###SBATCH --error=current_%A_%a.err
#SBATCH --time=800:00
##SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=owners

ml gsl
ml python/3.6.1
ml py-scipy/1.1.0_py36
ml viz
ml py-matplotlib/2.1.2_py36
ml cmake/3.13.1
ml py-numpy/1.17.2_py36
ml gcc/8.1.0

pwd
python3 -u sscpol/fit.py --method ps --blazar S5flare


