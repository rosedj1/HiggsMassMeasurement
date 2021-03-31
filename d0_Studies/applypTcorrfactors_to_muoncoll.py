"""MyMuonCollection pT Corrector

Open a pickled MyMuonCollection that has MyMuons in self.muon_ls.
Apply 

NOTE: Run this script in a dev session:
srun --partition=bigmem --mem=128gb --ntasks=1 --cpus-per-task=8 --time=08:00:00 --pty bash -i

Syntax: python this_script.py
Author: Jake Rosenzweig
Created: 2021-03-28
Updated: 2021-03-30
"""
import os
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, check_overwrite
from d0_Studies.d0_Utils.d0_fns import correct_muon_pT

use_GeoFit_algo = 1
corr_type = "GeoFit"  # AdHoc
overwrite = 0

inpkl_path_pTcorrfactors = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pTCorrFactors/GeoFit/GeoFitTcorrfact_derivedfromMC2016_3etabins_0p0eta2p4_newestformat.pkl"
inpkl_path_muoncoll      = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/MC2016DY_skim_fullstats_nogenmatching_0p01_d0_1000p0.pkl"
outpkl_path_muoncoll     = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/MC2016DY_skim_fullstats_nogenmatching_0p01_d0_1000p0_withGeoFitcorr.pkl"

assert corr_type in outpkl_path_muoncoll
if use_GeoFit_algo:
    assert corr_type in "GeoFit"
check_overwrite(outpkl_path_muoncoll, overwrite=overwrite)
muon_coll           = open_pkl(inpkl_path_muoncoll)
pT_corr_factor_dict = open_pkl(inpkl_path_pTcorrfactors)

muon_coll.apply_pTcorr_to_all_muons(pT_corr_factor_dict, use_GeoFit_algo=1, force_zero_intercept=True, verbose=True)
save_to_pkl(muon_coll, outpkl_path=outpkl_path_muoncoll, overwrite=overwrite)
# # Analyze each KB2D one at a time on SLURM.
# muon_coll.save_KB2Ds_separate_dcts(outdir=outdir_kb2d_dicts, file_prefix="", overwrite=overwrite, verbose=verbose)