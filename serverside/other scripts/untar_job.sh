#!/bin/bash
#SBATCH --job-name="ctd analyze"
#SBATCH --account=def-nereusvc
#SBATCH --mem-per-cpu=3900   #4 Gb
#SBATCH --time=06:00:00
#SBATCH --ntasks=4 # CPUs

tar -xf ~/projects/def-nereusvc/mdunphy/nemo_results/SalishSea1500/SalishSea1500-RUN141.tar -C ~/projects/def-nereusvc/goldford/nemo_results/SalishSea1500/

