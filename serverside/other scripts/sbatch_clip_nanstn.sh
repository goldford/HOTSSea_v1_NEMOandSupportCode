#!/bin/bash
#SBATCH --job-name="clip nanstn"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3000   #4 Gb = 3900 mb (max)
#SBATCH --time=01:00:00
#SBATCH --ntasks=10 # CPUs

ksh nanoose_clip_SS1500.ksh

