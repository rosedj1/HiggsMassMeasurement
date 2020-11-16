#!/bin/bash
#SBATCH --job-name=derive_corrections_0p0eta0p2
#SBATCH --output=job_%j.log
#SBATCH --error=job_%j_error.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rosedj1@ufl.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=80mb
#SBATCH --time=02:00:00
#SBATCH --account=avery
#SBATCH --qos=avery

pwd; hostname; date
source ~/.bash_profile
cd /blue/avery/rosedj1/HiggsMassMeasurement/
conda activate my_root_env
source setup.sh
echo "Packages loaded."

cd /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/
echo "Starting script..."
time python derive_pTcorrfactors_from_ggF_sample.py > /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Output/output_fullstats_0p0eta0p2.txt
echo "Script finished!"
