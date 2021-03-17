#!/bin/bash
#SBATCH --job-name=FILENAME
#SBATCH --output=OUTDIR/FILENAME.log
#SBATCH --error=/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/output/DeriveCorrections/2016DY/individualKB2Ditergaussfits/MC2016DY_kb2d_ETARANGE_PTRANGE_error.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rosedj1@ufl.edu
#SBATCH --time=04:00:00
#SBATCH --account=avery
#SBATCH --qos=avery-b
#SBATCH --mem=32gb
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
time python /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/submit_singlekb2d_itergaussfits.py
# time python derive_pTcorrfactors_from_Hmumu_sample.py > /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Output/DeriveCorrections/2018Hmumu/output_Hmumu_derivepTcorr_fullstats_synchwithXunwuetabins_50qd0reg.txt
echo "Script finished!"
