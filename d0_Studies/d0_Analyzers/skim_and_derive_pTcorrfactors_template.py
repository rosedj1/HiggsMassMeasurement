"""AdHoc pT Correction Deriver

Works for: ggH4mu, H2mu, DY

FIXME: This doc string is old!

Purpose:
    Extract muons from
    Put each muon into User-specified (eta, pT, q*d0) bins,
    Make dpT/pT and q*d0 distributions in each bin.
    Perform iterated Gaussian fits of dpT/pT distributions to describe sigma of core.
    Put all kinematic plots into a single PDF.
Syntax:
    python script.py > output.txt  (recommended)
    python script.py 
Notes:   
    Make sure to check all the parameters in "User Parameters".
    Should be used with Python 3.X.

    This code makes the following distributions, before/after pT corr:
        (1) muon pT dist.
        (2) m4mu dist. (DSCB fit?)
        (3) eta dist., for each eta bin
        (4) muon qd0 dist, per (eta, pT) bin
    The following plots are produced:
Author:  Jake Rosenzweig
Created: 2020-08-24
Updated: 2021-03-10
"""
import pickle
import os
import ROOT as r
import numpy as np
# Local imports.
from Utils_Python.Utils_Files import check_overwrite, make_dirs, save_to_pkl
from Utils_ROOT.Printer import CanvasPrinter
from ParticleCollections import MyMuonCollection
from Particles import MyMuon
from d0_Studies.kinematic_bins import equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum
from d0_Studies.d0_Utils.d0_cls import KinBin2D
from d0_Studies.d0_Utils.d0_fns import find_bin_edges_of_value
#----- User Parameters -----#
# Input.
infile_path = "/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/DY_JPsi_Upsilon/DY_2016.root"
prod_mode = "DY2mu"
# Output.
# outdir_pdf    = f"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/root_files/"
# outdir_pdf    = f"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/ggF_m4mu/"
outdir_pdf = "REPLACE_OUTDIR_PDF"
# outdir_pkl     = f"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/"  #MuonPickles/
outdir_pkl = "REPLACE_OUTDIR_PKL"
outfile_prefix = "REPLACE_JOB_NAME"

eta_ls = REPLACE_ETA_LS
# eta_ls = REPLACE_ETA_LS
# eta_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.20, 2.40]
# pT_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0]#, 1000.0] #bin_edges_pT_sevenfifths_to1000GeV_wholenum
# pT_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 1000.0] #bin_edges_pT_sevenfifths_to1000GeV_wholenum
# pT_ls = [5.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 1000.0] #bin_edges_pT_sevenfifths_to1000GeV_wholenum
# pT_ls = [5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 1000.0] # Merging 100-150-200-1000 GeV bins.
pT_ls = REPLACE_PT_LS

overwrite = REPLACE_OVERWRITE
verbose = REPLACE_VERBOSE
alarm_level = "warning" # "error"
show_plots_to_screen = False

make_dpTOverpT_vs_qd0_graphs = True
save_KinBin2D_pkl = True
save_pT_corr_pkl = True
# Decide to bin the qd0 axis or not.
do_unbinned_analysis = False
do_binned_analysis = True
binned_fit = REPLACE_binned_fit
use_data_in_xlim = REPLACE_use_data_in_xlim  # Only fit over data within x-window.
fit_whole_range_first_iter = REPLACE_fit_whole_range_first_iter  # False gives more consistent iterated fits.
min_muons_per_qd0_bin = REPLACE_min_muons_per_qd0_bin

max_n_evts = REPLACE_max_n_evts
print_out_every = REPLACE_print_out_every
regions = REPLACE_REGIONS       # Number of equal-entry regions along qd0 axis.
iters = REPLACE_ITERS
num_sigmas = 2.5  # Controls iterated fit range.
switch_to_binned_fit = REPLACE_switch_to_binned_fit

bins_dpTOverpT = 75
bins_qd0 = 100
x_lim_dpTOverpT = [-0.12, 0.12]  #[-0.08, 0.08]
x_lim_qd0 = [-0.01, 0.01]
#----- Functions -----#
def prep_files(eta_ls, pT_ls, 
               outfile_prefix, outdir_pdf, outdir_pkl, 
               regions, iters, overwrite=False):
    """Return the file paths for the produced files."""
    eta_min, eta_max, n_eta_bins = min(eta_ls), max(eta_ls), len(eta_ls)
    pT_min, pT_max, n_pT_bins    = min(pT_ls), max(pT_ls), len(pT_ls)
    bininfo  = f"{eta_min}eta{eta_max}w{n_eta_bins-1}bins"
    bininfo += f"_{pT_min:.0f}pT{pT_max:.0f}w{n_pT_bins-1}bins"
    new_file_name = f"{outfile_prefix}_{bininfo}_{regions}qd0reg_{min_muons_per_qd0_bin}minmu_{iters}itergausfit_{num_sigmas}sigs"
    new_file_name = new_file_name.replace('.', 'p')
    outpath_pdf = os.path.join(outdir_pdf, f"{new_file_name}.pdf")
    outpkl_path = os.path.join(outdir_pkl, f"{new_file_name}.pkl")
    outpkl_pT_corr_path = outpkl_path.replace(".pkl", "_pTcorrfactors.pkl")
    # See if these files already exist and overwrite if it's OK to do so.
    make_dirs(outdir_pdf, verbose)
    make_dirs(outdir_pkl, verbose)
    for f in [outpath_pdf, outpkl_path, outpkl_pT_corr_path]:
        check_overwrite(f, overwrite)
    return (outpath_pdf, outpkl_path, outpkl_pT_corr_path)

def main():
    assert not (do_binned_analysis and do_unbinned_analysis)
    if not show_plots_to_screen: 
        r.gROOT.SetBatch(True)
    # Prep your area.
    outpath_pdf, outpkl_path, outpkl_pT_corr_path = prep_files(eta_ls, pT_ls, outfile_prefix, outdir_pdf, outdir_pkl, 
                                                               regions, iters, overwrite=overwrite)
    printer = CanvasPrinter()
    printer.make_plots_pretty()
    # Begin analysis.
    mu_coll = MyMuonCollection(prod_mode=prod_mode)
    mu_coll.extract_muons(infile_path, prod_mode=prod_mode, n_evts=max_n_evts,
                                  print_out_every=print_out_every, eta_min=0.0, eta_max=2.4,
                                  pT_min=5, pT_max=200,
                                  d0_max=1, dR_max=0.002,
                                  do_mu_pT_corr=False,
                                  force_zero_intercept=False,
                                  pT_corr_factor_dict=None,
                                  use_GeoFit_algo=False,
                                  verbose=verbose)
    # Make inclusive plots.
    mu_coll.make_inclusive_kinematic_plots()
    # Sort muons into specified (eta, pT, qd0) bins.
    mu_coll.sort_muons(eta_ls, pT_ls, pT_corr_factor_dict=None,
                       n_bins_dpTOverpT=bins_dpTOverpT, x_lim_dpTOverpT=x_lim_dpTOverpT,
                       n_bins_qd0=bins_qd0, x_lim_qd0=x_lim_qd0, verbose=verbose)
    mu_coll.make_KinBin3Ds(regions=regions,
                           min_muons_per_qd0_bin=min_muons_per_qd0_bin, 
                           verbose=verbose)
    # Decide to do iter Gaus binning of qd0 or to do scatterplot of qd0.
    if do_binned_analysis:
        mu_coll.do_3D_iter_gaus_fits(binned_fit=binned_fit, bins_dpTOverpT=bins_dpTOverpT, bins_qd0=bins_qd0, 
                                     x_lim_dpTOverpT=x_lim_dpTOverpT, x_lim_qd0=x_lim_qd0, 
                                     fit_whole_range_first_iter=fit_whole_range_first_iter,
                                     iters=iters,
                                     num_sigmas=num_sigmas,
                                     switch_to_binned_fit=switch_to_binned_fit,
                                     alarm_level=alarm_level,
                                     verbose=verbose, use_data_in_xlim=use_data_in_xlim)
        if make_dpTOverpT_vs_qd0_graphs:
            mu_coll.make_KinBin2D_graphs()
            mu_coll.make_all_multigraphs(eta_ls)  # Appends to self.multigraph_ls
            # Get dpT/pT vs. qd0 graphs.
            # mu_coll.collect_KinBin2D_plots()
        # Get dpT/pT IterGaussFit plots and qd0 dists.
        # mu_coll.collect_KinBin3D_plots()

        # printer.make_pdf_of_plots(mu_coll.get_all_plots(), outpath_pdf)

    elif do_unbinned_analysis:
        # Make a scatterplot of dpT/pT vs. qd0 values for each muon. 
        # Make all plots.
        raise RuntimeError("This section needs to be updated.")
        # mu_coll.make_inclusive_kinematic_plots()
        mu_coll.fill_KinBin_hists()
        mu_coll.make_KinBin2D_graphs()
        mu_coll.collect_all_plots()
        # Print plots into single PDF.
        print(f"Drawing all plots to:\n{outpath_pdf}")
        printer.make_pdf_of_plots(mu_coll.plot_ls, outpath_pdf)
        # mu_coll.save_to_pkl(mu_coll.pT_corr_factor_dict, outpkl_pT_corr_path)

    # Draw all plots.
    printer.canv.Print(outpath_pdf + "[")
    mu_coll.draw_all_multigraphs(outpath_pdf, printer)
    # "graphs", "dpT/pT hists", "qd0 hists"
    all_plots = mu_coll.get_all_kb3d_plots("dpT/pT iterfit", "qd0 hists") + mu_coll.get_all_kb2d_plots("dpT/pT hists", "qd0 hists") + mu_coll.hist_inclusive_ls
    printer.draw_plots(all_plots, outpath_pdf)
    # printer.draw_plots(mu_coll.hist_inclusive_ls + mu_coll.kinbin2d_plot_ls + mu_coll.kinbin3d_plot_ls,
    #                    outpath_pdf)
    printer.canv.Print(outpath_pdf + "]")

    if save_pT_corr_pkl:
        mu_coll.make_pT_corr_dict()
        save_to_pkl(mu_coll.pT_corr_factor_dict, outpkl_pT_corr_path)
    if save_KinBin2D_pkl:
        mu_coll.overwrite_kb2d_muon_info()
        # mu_coll.overwrite_kb3d_muon_info()
        # save_to_pkl(mu_coll.KinBin2D_dict, outpkl_path)
        save_to_pkl(mu_coll, outpkl_path)

if __name__ == "__main__":
    main()