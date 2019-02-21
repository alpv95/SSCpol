#!/bin/bash
#
#SBATCH --job-name=testALP
#SBATCH --time=7:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=iric
#SBATCH --nodes=1
#SBATCH --output=current.out

module load gcc
module load python/3.6.1
module load py-numpy/1.14.3_py36
module load py-pytorch/1.0.0_py36
module load viz
module load py-matplotlib/2.1.2_py36

srun python3 automateSED.py
