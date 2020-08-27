#!/bin/bash
#SBATCH --job-name=jupyter_python2
#SBATCH --output=notebook_jupyter_python2_%j.log
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=08:00:00
date;hostname;pwd
 
module add jupyter

cd /blue/avery/rosedj1/HiggsMassMeasurement/CMSSW_10_2_18/src/
eval `scramv1 runtime -sh`
cd /blue/avery/rosedj1/

source /blue/avery/rosedj1/HiggsMassMeasurement/setup.sh

launch_jupyter_notebook
