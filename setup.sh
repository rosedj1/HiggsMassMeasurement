#!/bin/bash

export PYTHONPATH=${PYTHONPATH}:${PWD}
export PYTHONPATH=${PYTHONPATH}:${PWD}/d0_Studies
export PYTHONPATH=$PYTHONPATH:/blue/avery/rosedj1/HiggsMassMeasurement
export PYTHONPATH=$PYTHONPATH:/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies
#export PYTHONPATH=$PYTHONPATH:${PWD}/leptonPtErrorCorrector/
#export PYTHONPATH=$PYTHONPATH:${PWD}/leptonPtErrorCorrector/doCorrection/
#export PYTHONPATH=$PYTHONPATH:${PWD}/leptonPtErrorCorrector/doCorrection/tmp/
#export PYTHONPATH=$PYTHONPATH:${PWD}/doClosure/
#export PYTHONPATH=$PYTHONPATH:${PWD}/doClosure/ZClosure/
#export PYTHONPATH=$PYTHONPATH:/home/rosedj1/HiggsMeasurement/CMSSW_10_2_15/src/HiggsMassMeasurement/
#export PYTHONPATH=$PYTHONPATH:/home/rosedj1/HiggsMeasurement/CMSSW_10_2_15/src/HiggsMassMeasurement/d0_Studies/

#export MPLCONFIGDIR=/home/rosedj1/.config/matplotlib:$MPLCONFIGDIR
export MPLCONFIGDIR=/home/rosedj1/.config/matplotlib/mpl_configdir/stylelib/:$MPLCONFIGDIR

export PATH=${PATH}:${PWD}/bin
export BASE_PATH=${PWD}

#cd /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_15/src/
#eval `scramv1 runtime -sh`
#cd -
