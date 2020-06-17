#!/bin/bash
#SBATCH --job-name=jupyter_python3
#SBATCH --output=notebook_jupyter_python3_%j.log
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --time=08:00:00
date;hostname;pwd
 
module add jupyter

cd /ufrc/avery/rosedj1/HiggsMassMeasurement/CMSSW_10_2_18/src/
eval `scramv1 runtime -sh`
cd /ufrc/avery/rosedj1/HiggsMassMeasurement/

ml python3
ml texlive  # This lets us have LaTeX font in Matplotlib plots.
#ml conda
#conda init bash
#conda activate my_root_env

source setup.sh

launch_jupyter_notebook
