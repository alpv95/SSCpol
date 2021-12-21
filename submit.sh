#!/bin/bash
#
#SBATCH --job-name=S5
#SBATCH --output=slurm/%A.out
####SBATCH --error=current_%A_%a.err
#SBATCH --time=800:00
##SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=owners

python3 -u fit.py --method standard


