#!/bin/bash
#SBATCH --job-name=REPLACE_JOB_NAME_REPLACE_ETA_NAME
#SBATCH --output=REPLACE_OUTDIR_TXT/REPLACE_JOB_NAME_REPLACE_ETA_NAME.log
#SBATCH --error=REPLACE_OUTDIR_TXT/REPLACE_JOB_NAME_REPLACE_ETA_NAME_error.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rosedj1@ufl.edu
#SBATCH --time=12:00:00
#SBATCH --account=avery
#SBATCH --qos=avery
#SBATCH --mem=4gb
#SBATCH --partition=hpg2-compute
#SBATCH --nodes=2

pwd; hostname; date
source ~/.bash_profile
cd /blue/avery/rosedj1/HiggsMassMeasurement/
conda activate my_root_env
source setup.sh
echo "Packages loaded."

cd /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/
echo "Starting script..."
time python REPLACE_NEW_SCRIPT
echo "Script finished!"

