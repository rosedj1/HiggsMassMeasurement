"""Template for Iterative Gaussian Fits on Rochester Corr. vs. Non-RC samples

Requires an input pickled dict of KinBin2Ds.
- You can make the input dict from roch_vs_noroch_kb2dmaker_inclusivehistplotter.py

This is a template file whose values (REPLACE_x) get replaced by a SLURM
submission script.

The iterative fits are computationally intensive after about 3 fits,
so you should use the following SLURM submission file to run this one:
- d0_Studies/d0_Analyzers/submit_to_slurm_inbatch.py

Author: Jake Rosenzweig
Updated: 2021-03-04
"""
from Utils_Python.Utils_Files import open_pkl, make_dirs
from Utils_Python.Utils_Files import check_overwrite
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid
from Utils_ROOT.ROOT_StatsAndFits import RooFit_iterative_gaus_fit
from d0_Studies.RochCorrAnalyzers.roch_vs_noroch_itergausfit_local import run_script, prep_area
import ROOT
import os
import gc
#----- User Switches -----#
overwrite = REPLACE_OVERWRITE
iters = REPLACE_ITERS
verbose = REPLACE_VERBOSE
fit_whole_range_first_iter = REPLACE_fit_whole_range_first_iter  # False gives more consistent fits (with no outlier data).
use_data_in_xlim = REPLACE_use_data_in_xlim
binned_fit = REPLACE_binned_fit
switch_to_binned_fit = REPLACE_switch_to_binned_fit
# eta_ls = equal_entry_bin_edges_eta_mod1_wholenum #[0.0, 0.9, 1.7, 2.4]
# pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum#[5, 25, 50, 100]
inpath_pkl = "REPLACE_INPATH_PKL"
# inpath_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/pickles/test/fixitergaussfits_test01_3kbds.pkl"
# Outpath info is controlled by d0_Studies/d0_Analyzers/submit_to_slurm_inbatch.py
filename   = "REPLACE_JOB_NAME"
outdir_pkl = "REPLACE_OUTDIR_PKL"
outdir_txt = "REPLACE_OUTDIR_TXT"  # os.path.join("REPLACE_OUTDIR_TXT", filename)
outdir_pdf = "REPLACE_OUTDIR_PDF"  # os.path.join("REPLACE_OUTDIR_PDF", filename)
eta_range = REPLACE_ETA_LS  # [0.0, 0.2]

if __name__ == "__main__":
    outpath_pkl, outpath_txt, outpath_pdf = prep_area(eta_range, filename, iters, outdir_pkl, outdir_txt, outdir_pdf, overwrite)
    run_script(inpath_pkl, eta_range, binned_fit, switch_to_binned_fit, iters,
               fit_whole_range_first_iter, use_data_in_xlim, outpath_pdf, outpath_pkl)