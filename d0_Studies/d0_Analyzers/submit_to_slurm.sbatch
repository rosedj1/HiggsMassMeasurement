#!/bin/bash
#SBATCH --job-name=synch50qd0
#SBATCH --output=job_%j.log
#SBATCH --error=job_%j_error.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rosedj1@ufl.edu
#SBATCH --time=12:00:00
#SBATCH --account=avery
#SBATCH --qos=avery
#SBATCH --mem=16gb
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
time python derive_pTcorrfactors_from_Hmumu_sample.py > /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Output/DeriveCorrections/2018Hmumu/output_Hmumu_derivepTcorr_fullstats_synchwithXunwuetabins_50qd0reg.txt
echo "Script finished!"
