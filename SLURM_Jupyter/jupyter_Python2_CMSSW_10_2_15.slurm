#!/bin/bash
#SBATCH --job-name=jupyter_python2
#SBATCH --output=notebook_jupyter_python2_%j.log
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --time=06:00:00
date;hostname;pwd
 
module add jupyter

cd /ufrc/avery/rosedj1/HiggsMassMeasurement/CMSSW_10_2_15/src/
eval `scramv1 runtime -sh`
cd /ufrc/avery/rosedj1/

launch_jupyter_notebook
