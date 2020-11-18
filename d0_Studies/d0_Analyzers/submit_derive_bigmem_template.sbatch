#!/bin/bash
#SBATCH --job-name=2017derive_fullstats_REPLACE_ETA_NAME
#SBATCH --output=job_%j.log
#SBATCH --error=job_%j_error.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rosedj1@ufl.edu
#SBATCH --time=12:00:00
#SBATCH --account=avery
#SBATCH --qos=avery
#SBATCH --mem=2gb
#SBATCH --partition=bigmem
#SBATCH --nodes=1

pwd; hostname; date
source ~/.bash_profile
cd /blue/avery/rosedj1/HiggsMassMeasurement/
conda activate my_root_env
source setup.sh
echo "Packages loaded."

cd /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Plotters/
echo "Starting script..."
# time python plot_m4mu_kinematics_withpTcorr.py > /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Output/ApplyCorrections/output_dpTOverpT_beforeafterpTcorr_fullstats.txt
time python REPLACE_NEW_FILE > /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Output/DeriveCorrections/2017ggH/output_fullstats_REPLACE_ETA_NAME.txt
echo "Script finished!"
