#!/bin/bash
#SBATCH --job icon_system2
# send email about job start and end
#SBATCH --mail-user=bka28@pitt.edu
#SBATCH --mail-type=ALL

#SBATCH --nodes=1
#SBATCH --partition=dept_gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint="L40"


#run the MD in the conda environment "struct"
conda run -n md python dcdTools.py







