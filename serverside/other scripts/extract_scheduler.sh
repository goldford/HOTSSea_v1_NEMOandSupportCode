#!/bin/bash
#SBATCH --job-name="extract custom GO"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3900   #4 Gb
#SBATCH --time=08:00:00 # 
#SBATCH --ntasks=6 # CPUs

# Load necessary modules
#module load python/3.8.10

# Activate python environment
#source activate your_env

# Run your Python script
#python extract_results_1_toweekly.py
#python extract_results_2_anomclim.py
python extract_results_3_anomtrends.py
#python extract_results_4_rawtrends.py