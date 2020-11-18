#!/bin/bash
#SBATCH --job-name=fullstats
#SBATCH --output=job_%j.log
#SBATCH --error=job_%j_error.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rosedj1@ufl.edu
#SBATCH --time=12:00:00
#SBATCH --account=avery
#SBATCH --qos=avery
#SBATCH --mem=32gb
#SBATCH --partition=bigmem
#SBATCH --nodes=2

pwd; hostname; date
source ~/.bash_profile
cd /blue/avery/rosedj1/HiggsMassMeasurement/
conda activate my_root_env
source setup.sh
echo "Packages loaded."

cd /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Plotters/
echo "Starting script..."
time python plot_m4mu_kinematics_withpTcorr.py > /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Output/ApplyCorrections/output_dpTOverpT_beforeafterpTcorr_fullstats.txt
echo "Script finished!"
