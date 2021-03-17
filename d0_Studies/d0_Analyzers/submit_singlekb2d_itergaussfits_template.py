"""Submit Iterated Gaussian Fit to SLURM

Purpose: This script opens a pickled KinBin2D with MyMuons
sorted into all the KB3Ds.
Iterated Gaussian fits (IGFs) are performed on the dpT/pT dist of each KB3D.
The fit results are saved to the KB3D.
The processed KB2D is then pickled into a new file.

Syntax: python this_script.py
Author: Jake Rosenzweig
Created: 2021-03-14  # Happy pi day!
Updated: 2021-03-15
"""
import os
from Utils_ROOT.Printer import CanvasPrinter
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, check_overwrite, make_dirs
from d0_Studies.d0_Utils.d0_fns import print_header_message
from d0_Studies.d0_Analyzers.slurm_inbatch_derive_pTcorrfactors import make_name_from_ls
from pprint import pprint

overwrite = 1

bins_dpTOverpT = 100
bins_qd0 = 100
x_lim_dpTOverpT = [-0.4, 0.4]
x_lim_qd0 = [-0.01, 0.01] 
binned_fit = False
fit_whole_range_first_iter = False 
iters = 5
num_sigmas = 2.5
switch_to_binned_fit = 999999
verbose = True
alarm_level = "warning"
use_data_in_xlim = True

eta_range = REPLACE_ETA_LS #[0.4, 0.6]
pT_range = REPLACE_PT_RANGE#[20.0, 27.0]

# inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/kb2d_dictsnofitinfo/MC2016DY_fullstats_muoncoll_withkb3dbins__0p0eta0p2*.pkl"
inpkl_path = "REPLACE_INPKL_PATH"
outpkl_path = "REPLACE_OUTPATH_PKL"
# outpkl_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/kb2d_dictswithfitinfo/test/REPLACE_JOB_NAME_{eta_name}_{pT_name}.pkl"

def loop_over_kb3d(kb2d):
    """Do IGFs for all KB3Ds in this KB2D."""
    for kb3d in kb2d.KinBin3D_dict.values():
        kb3d.analyze_KinBin3D(bins_dpTOverpT=bins_dpTOverpT, bins_qd0=bins_qd0,
                            x_lim_dpTOverpT=x_lim_dpTOverpT, x_lim_qd0=x_lim_qd0, 
                            binned_fit=binned_fit, fit_whole_range_first_iter=fit_whole_range_first_iter, 
                            iters=iters, num_sigmas=num_sigmas, 
                            switch_to_binned_fit=switch_to_binned_fit,  
                            verbose=verbose, alarm_level=alarm_level, 
                            use_data_in_xlim=use_data_in_xlim)
        print_header_message(f"KB3D done: eta={kb3d.eta_range} pT={kb3d.pT_range}")

def prep_area(outpkl_path, overwrite=False):
    """Make dirs, check overwrite."""
    outdir = os.path.dirname(outpkl_path)
    make_dirs(outdir)
    check_overwrite(outpkl_path, overwrite)

def main():
    printer = CanvasPrinter()
    printer.make_plots_pretty()
    prep_area(outpkl_path, overwrite=overwrite)
    kb2d = open_pkl(inpkl_path)
    print(f"Opened pkl:\n{inpkl_path}")
    kb2d.is_qd0_binned = True
    kb2d.overwrite_muon_info(delete_all=True)  # Muons now stored in KB3D.
    loop_over_kb3d(kb2d)
    # Make dpT/pT and dpT/pT * 1/<pT> vs. qd0 plots.
    kb2d.make_dpTOverpT_graph(color=4, do_fit=True, scale_by_1divpT=False)
    kb2d.make_dpTOverpT_graph(color=4, do_fit=True, scale_by_1divpT=True)
    for kb3d in kb2d.KinBin3D_dict.values():
        kb3d.overwrite_muon_info(delete_all=True)  # Don't need muons anymore.
    save_to_pkl(kb2d, outpkl_path)

if __name__ == "__main__":
    main()