"""SLURM Submitter for MyMuon Reco pT Corrector

Apply pT corrections stored in a '.pkl' to a list of individual pickled KB2Ds.
The muons will the store their pT_reco and pT_corr values.

NOTE: This script doesn't yet produce m4mu_corr.
The only way to evaluate m2mu_corr or m4mu_corr is by extracting muons from a
file. Thus far you can't open up a pickled MyMuonCollection and apply pT corr
to uncorrected m4mu vals.
Idea: class Event, which can store all 4 muons per event (hence m4mu).

This script controls:
d0_Studies/d0_Analyzers/singlekb2ditergaussfit_slurm_template.py

Syntax: python this_script.py
Author: Jake Rosenzweig
Created: <2021-04-12
Updated: 2021-04-12
"""
import os
import shutil
from glob import glob
from natsort import natsorted
from Utils_Python.Utils_Files import open_pkl, check_overwrite, replace_value, make_dirs
from Utils_Python.SlurmManager import SLURMSubmitter
from d0_Studies.kinematic_bins import equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum

verbose = 1
overwrite = 0
delete_kb2d_muon_ls = 1
corrections_derived_from = "DY"
method = "AdHoc"

inpkl_path_pTcorrfactordict = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/muoncoll_itergaussfitsonKB3Ds_redo_pTcorrfactors.pkl"
# inpkl_path_pTcorrfactordict = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/muoncollwithfitstats/muoncoll_itergaussfitsonKB3Ds_pTcorrfactors.pkl"
inpkl_kb2d_path_glob = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/skim_fullstats_verify/pickles/kb2d_dicts/*.pkl"
# inpkl_kb2d_path_glob = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/skim_fullstats_verify/pickles/kb2d_dicts/*.pkl"
# inpkl_kb2d_path_glob = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p01_d0_1000p0/*eta*_*pT*0.pkl"
fullpath_main_script = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/singlekb2ditergaussfit_slurm_template.py"

new_filename_prefix = "" #"unbinnedfit_widerwindow_fitwholerangefirstiter"
outdir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/ApplyCorr/MC2016DY2mu/individKB2Ds_beforeaftercorr/"
# outcopies_dir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/beforeafterGeoFitCorr/0p01_d0_1000p0/copies/"
# outtxt_dir    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/beforeafterGeoFitCorr/0p01_d0_1000p0/output/"

# eta_binedge_ls = equal_entry_bin_edges_eta_mod1_wholenum
# pT_binedge_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum

bins_dpTOverpT = 200
bins_qd0 = 200
x_lim_dpTOverpT = [-0.5, 0.5]
x_lim_qd0 = [-0.01, 0.01]
use_smart_window = 1  # Disregards user-specified x limits.
iters = 5

switch_to_binned_fit = 99E9

# SLURM directives.
partition = "hpg2-compute" #"bigmem"
mem = (24, 'gb')
nodes = 1
burst = False
time = "01:00:00"  # hh:mm:ss

#--- Main ---#
outcopies_dir = os.path.join(outdir, "copies/")
outtxt_dir    = os.path.join(outdir, "output/")
outpkl_dir    = os.path.join(outdir, "pickles/")
assert corrections_derived_from in inpkl_path_pTcorrfactordict
assert "." not in new_filename_prefix
# Perform iterated Gaussian fits on each KinBin2D.
kb2d_pkl_ls = glob(inpkl_kb2d_path_glob)
kb2d_pkl_ls = natsorted(kb2d_pkl_ls)

assert len(kb2d_pkl_ls) > 0, "[ERROR] No files found."
# template_basename = template.rstrip(".py")
for inpkl_path in kb2d_pkl_ls:
    # Create main copy.
    base = os.path.basename(fullpath_main_script).rstrip(".py")
    kb2d_bin_key = os.path.basename(inpkl_path).rstrip(".pkl")
    if len(new_filename_prefix) > 0:
        base = f"{new_filename_prefix}_{base}"
    file_name_copy = f"{base}_copy_{kb2d_bin_key}"
    template_copy = os.path.join(outcopies_dir, f"{file_name_copy}.py")
    make_dirs(outcopies_dir)
    make_dirs(outtxt_dir)
    make_dirs(outpkl_dir)
    shutil.copyfile(fullpath_main_script, template_copy)
    outpkl_suffix = f"{kb2d_bin_key}_withkb2dfits.pkl"
    outpkl_filename = f"{new_filename_prefix}_{outpkl_suffix}" if len(new_filename_prefix) > 0 else outpkl_suffix
    outpkl_path = os.path.join(outpkl_dir, outpkl_filename)
    # Modify main copy.
    replace_value("OVERWRITE", overwrite, template_copy)
    replace_value("VERBOSE", verbose, template_copy)
    replace_value("INPKL_PATH", inpkl_path, template_copy)
    replace_value("INPKL_PTCORRFACTORDICT", inpkl_path_pTcorrfactordict, template_copy)
    replace_value("OUTPKL_PATH", outpkl_path, template_copy)
    replace_value("DELETE_KB2D_MUON_LS", delete_kb2d_muon_ls, template_copy)
    replace_value("SWITCH2BINNED", switch_to_binned_fit, template_copy)
    # replace_value("ETA_BINEDGE_LS", eta_binedge_ls, template_copy)
    # replace_value("PT_BINEDGE_LS", pT_binedge_ls, template_copy)
    replace_value("BINS_DPTOVERPT", bins_dpTOverpT, template_copy)
    replace_value("BINS_QD0", bins_qd0, template_copy)
    replace_value("X_LIM_DPTOVERPT", x_lim_dpTOverpT, template_copy)
    replace_value("X_LIM_QD0", x_lim_qd0, template_copy)
    replace_value("ITERS", iters, template_copy)
    replace_value("USE_SMART_WINDOW", use_smart_window, template_copy)
    
    sbmtr = SLURMSubmitter(verbose=verbose)
    sbmtr.prep_directives(
        job_name=file_name_copy,
        output_txt=os.path.join(outtxt_dir, f"{file_name_copy}.log"),
        email="rosedj1@ufl.edu",
        time=time,
        acct="avery",
        burst=burst,
        mem=mem,
        partition=partition,
        nodes=nodes,
    )
    cmdtup = (f"python {template_copy}")
    fullpath_slurm_copy = os.path.join(outcopies_dir, f"{file_name_copy}.sbatch")
    result = sbmtr.make_slurm_script(fullpath_slurm_copy, cmdtup, overwrite=overwrite)
    if result == 0:
        sbmtr.submit_script(fullpath_slurm_copy)
    else:
        raise RuntimeError("SLURM script not submitted.")

print(f"[INFO] Copies stored at:\n{outcopies_dir}")
print(f"[INFO] Output stored at:\n{outtxt_dir}")
print(f"[INFO] Pickle stored at:\n{outpkl_dir}")