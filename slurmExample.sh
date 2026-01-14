#!/bin/bash
#SBATCH --job <jobName>
#SBATCH --mail-user=<email>
#SBATCH --mail-type=ALL

#SBATCH --nodes=1
#SBATCH --partition=dept_gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint="L40"

#run the MD in the conda environment "md"
conda run -n md python dcdTools.py
