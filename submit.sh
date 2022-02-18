#!/bin/bash
#
#SBATCH --job-name=BlazarFit
#SBATCH --output=slurm/%A.out
#SBATCH --time=800:00
##SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=31
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=owners

python3 -u fit.py --blazar S5flare --method cross_entropy --nblocks 37 --nprocs 30 --rand_gamma 0
