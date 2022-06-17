"""Submit Iterated Gaussian Fits to SLURM

TODO:
[ ] Get rid of many User Parameters and put them into:
HiggsMassMeasurement/d0_Studies/d0_Analyzers/slurm_inbatch_singlekb2d_itergaussfitsonkb3ds.py

Purpose:
This script opens a pickled KinBin2D with MyMuons in self.muon_ls.
Muons are sorted into KB3Ds.
Iterated Gaussian fits (IGFs) are performed on the dpT/pT dist of each KB3D.
The fit results are saved to the KB3D.
The processed KB2D is then pickled into a new file.

Syntax: python this_script.py
Author: Jake Rosenzweig
Created: 2021-03-14  # Happy pi day!
Updated: 2021-04-14
"""
import os
from Utils_ROOT.Printer import CanvasPrinter
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, check_overwrite, make_dirs
from Utils_Python.printing import announce
from d0_Studies.d0_Analyzers.slurm_inbatch_derive_pTcorrfactors import make_name_from_ls
from pprint import pprint

verbose = VERBOSE
overwrite = OVERWRITE

bins_dpTOverpT = 100
x_lim_dpTOverpT = [-0.4, 0.4] #None

bins_qd0 = 100
x_lim_qd0 = [-0.01, 0.01]

bins_1OverpT = 100
x_lim_1OverpT = [0, 0.2] #None

binned_fit = BINNED_FIT
fit_whole_range_first_iter = False
use_smart_window = USE_SMART_WINDOW
iters = ITERS
num_sigmas = 2.5
use_data_in_xlim = True
# auto_range = True
switch_to_binned_fit = 999999999
alarm_level = "warning"
fit_with_zero_interc = FIT_WITH_ZERO_INTERC

regions = REGIONS  # Number of q*d0 bins per KB2D.
min_muons_per_qd0_bin = MIN_MUONS_PER_QD0_BIN

delete_kb3d_muons = True

eta_range = REPLACE_ETA_LS #[0.4, 0.6]
pT_range = REPLACE_PT_RANGE#[20.0, 27.0]

# inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/kb2d_dictsnofitinfo/MC2016DY_fullstats_muoncoll_withkb3dbins__0p0eta0p2*.pkl"
inpkl_path = "REPLACE_INPKL_PATH"
outpkl_path = "REPLACE_OUTPATH_PKL"
# outpkl_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/kb2d_dictswithfitinfo/test/REPLACE_JOB_NAME_{eta_name}_{pT_name}.pkl"

def do_itergaussfits_on_kb3ds(kb2d):
    """Do IGFs for all KB3Ds in this KB2D."""
    for kb3d in kb2d.KinBin3D_dict.values():
        # Do dpT/pT and 1/pT iterated Gaussian fits.
        kb3d.analyze_KinBin3D(bins_dpTOverpT=bins_dpTOverpT, bins_qd0=bins_qd0,
                            x_lim_dpTOverpT=x_lim_dpTOverpT, x_lim_qd0=x_lim_qd0,
                            bins_1OverpT=bins_1OverpT, x_lim_1OverpT=x_lim_1OverpT,
                            binned_fit=binned_fit, fit_whole_range_first_iter=fit_whole_range_first_iter, 
                            iters=iters, num_sigmas=num_sigmas, 
                            switch_to_binned_fit=switch_to_binned_fit,  
                            verbose=verbose, alarm_level=alarm_level, 
                            use_data_in_xlim=use_data_in_xlim,
                            use_smart_window=use_smart_window)
        
        announce(f"KB3D done: eta={kb3d.eta_range} pT={kb3d.pT_range}")

def prep_area(outpkl_path, overwrite=False):
    """Make dirs, check overwrite."""
    outdir = os.path.dirname(outpkl_path)
    make_dirs(outdir)
    check_overwrite(outpkl_path, overwrite)

def main():
    prep_area(outpkl_path, overwrite=overwrite)
    # Open MyMuonCollection with MyMuons separated into KB2Ds.
    # Then save each KB2D as its own pickled. dict.
    print(f"...Opening pickled KB2D:\n  {inpkl_path}")
    kb2d = open_pkl(inpkl_path)
    kb2d.is_qd0_binned = True
    # Make KB3Ds for this KB2D.
    print(f"...Making equal-entry KB3Ds for KB2D: eta={kb2d.eta_range}, pT={kb2d.pT_range}")
    kb2d.make_empty_equalentry_KinBin3Ds(regions, algo=("at_least", min_muons_per_qd0_bin), 
                                         verbose=verbose, title_friendly=False)
    print("...Placing muons into KB3Ds.")
    kb2d.store_muon_info_in_KinBin3Ds(title_friendly=False, verbose=verbose)
    kb2d.overwrite_muon_info(delete_all=True)  # Muons now stored in KB3D.
    # 
    printer = CanvasPrinter()
    printer.make_plots_pretty()
    do_itergaussfits_on_kb3ds(kb2d)
    # Make dpT/pT and dpT/pT * 1/<pT> vs. qd0 plots.
    kb2d.make_dpTOverpT_graph(color=4, do_fit=True, scale_by_1divpT=False, fit_with_zero_interc=fit_with_zero_interc)
    kb2d.make_dpTOverpT_graph(color=4, do_fit=True, scale_by_1divpT=True, fit_with_zero_interc=fit_with_zero_interc)
    kb2d.make_dpTOverpT_graph(color=4, do_fit=True, scale_by_avgOf1divpT=True, fit_with_zero_interc=fit_with_zero_interc)
    # kb2d.make_dpTOverpT_graph(color=4, do_fit=True, scale_by_muOf1divpT=True, fit_with_zero_interc=fit_with_zero_interc)
    # Also make dpT/pT * <1/pT> vs. qd0 plots (notice where < > is!).
    if delete_kb3d_muons:
        for kb3d in kb2d.KinBin3D_dict.values():
            kb3d.overwrite_muon_info(delete_all=True)  # Don't need muons anymore.
    save_to_pkl(kb2d, outpkl_path, overwrite=overwrite)

if __name__ == "__main__":
    main()