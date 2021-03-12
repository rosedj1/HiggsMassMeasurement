#!/bin/bash
#SBATCH --job-name=MC2017DY_muons
#SBATCH --output=/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/output/DeriveCorrections/2017DY/MC2017DY_muons.log
#SBATCH --error=/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/output/DeriveCorrections/2017DY/MC2017DY_muons_error.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rosedj1@ufl.edu
#SBATCH --time=8:00:00
#SBATCH --account=avery
#---SBATCH --qos=avery-b
#SBATCH --qos=avery
#SBATCH --mem=128gb
#SBATCH --partition=bigmem
#SBATCH --nodes=4

pwd; hostname; date
source ~/.bash_profile
cd /blue/avery/rosedj1/HiggsMassMeasurement/
conda activate my_root_env
source setup.sh
echo "Packages loaded."

cd /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/
echo "Starting script..."
time python -m cProfile -o output_2017.prof skim_sample.py
echo "Script finished!"