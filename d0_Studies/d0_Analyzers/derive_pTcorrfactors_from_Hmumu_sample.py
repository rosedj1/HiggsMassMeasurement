"""Ad Hoc pT Correction Analyzer

FIXME: This doc string is old!

Purpose:
    Select "good" H->ZZ->4mu events.
    Put each muon into User-specified (eta, pT, q*d0) bins,
    Make dpT/pT and q*d0 distributions in each bin.
    Perform iterated Gaussian fits of dpT/pT distributions to describe sigma of core. 
    all kinematic plots into a single PDF.
Syntax:
    python script.py > output.txt  (recommended)
    python script.py 
Notes:
    This code runs on a H->2mu sample.
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
Updated: 2020-12-16
"""
import pickle
import os
import ROOT as r
import numpy as np
# Local imports.
from Utils_Python.Utils_Files import check_overwrite, makeDirs
from Utils_ROOT.Printer import CanvasPrinter
from ParticleCollections import MyMuonCollection
from Particles import MyMuon
from d0_Studies.kinematic_bins import (equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum)
from d0_Studies.d0_Utils.d0_cls import KinBin2D
from d0_Studies.d0_Utils.d0_fns import find_bin_edges_of_value
#----- User Parameters -----#
# Input.
infile_path = f"/cmsuf/data/store/user/t2/users/rosedj1/Higgs/HiggsMassMeas/NTuples/Hmumu/MC2018Hmumu_liteskim_fullstats.root"
# Output.
# outpdf_dir    = f"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/root_files/"
# outpdf_dir    = f"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/ggF_m4mu/"
outpdf_dir = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Plots/derivepTcorrplots/Hmumu_2018/"
# outpkl_dir     = f"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/"  #MuonPickles/
outpkl_dir = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Pickles/2018/Hmumu/"
# Don't add any extension, like '.pdf' or '.pkl':
outfile_prefix = f"MC2018_Hmumu_derivepTcorr_fullstats_synchwithXunwuetabins"

# eta_ls = equal_entry_bin_edges_eta_mod1_wholenum[0:2]
eta_ls = [0.0, 0.9, 1.7, 2.4]
# eta_ls = REPLACE_ETA_LS
# eta_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.20, 2.40]
# pT_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0]#, 1000.0] #bin_edges_pT_sevenfifths_to1000GeV_wholenum
# pT_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 1000.0] #bin_edges_pT_sevenfifths_to1000GeV_wholenum
# pT_ls = [5.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 1000.0] #bin_edges_pT_sevenfifths_to1000GeV_wholenum
# pT_ls = [5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 1000.0] # Merging 100-150-200-1000 GeV bins.

pT_ls = [20.0, 1000.0]  # Must be floats.

overwrite = False
verbose = True
alarm_level = "warning" # "error"
show_plots_to_screen = False

make_dpTOverpT_vs_qd0_graphs = True
save_KinBin2D_pkl = True
save_pT_corr_pkl = True
# Decide to bin the qd0 axis or not.
# Choose one or the other.
do_unbinned_analysis = False
do_binned_analysis = True
fit_whole_range_first_iter = False
min_muons_per_qd0_bin = 100

max_n_evts = -1
print_out_every = 100000
regions = 50       # Number of equal-entry regions along qd0 axis.
iters = 5
num_sigmas = 2.5  # Controls iterated fit range.
switch_to_binned_fit = 5000

bins_dpTOverpT = 75
bins_qd0 = 100
# x_lim_dpTOverpT = [-0.1, 0.1]
x_lim_dpTOverpT = [-0.10, 0.10]
x_lim_qd0 = [-0.01, 0.01]
#----- Functions -----#
def prep_files(eta_ls, pT_ls, 
               outfile_prefix, outpdf_dir, outpkl_dir, 
               regions, iters, overwrite=False):
    """Return the file paths for the produced files."""
    eta_min, eta_max, n_eta_bins = min(eta_ls), max(eta_ls), len(eta_ls)
    pT_min, pT_max, n_pT_bins    = min(pT_ls), max(pT_ls), len(pT_ls)
    bininfo  = f"{eta_min}eta{eta_max}w{n_eta_bins-1}bins"
    bininfo += f"_{pT_min:.0f}pT{pT_max:.0f}w{n_pT_bins-1}bins"
    new_file_name = f"{outfile_prefix}_{bininfo}_{regions}qd0reg_{min_muons_per_qd0_bin}minmu_{iters}itergausfit_{num_sigmas}sigs"
    new_file_name = new_file_name.replace('.', 'p')
    outpdf_path = os.path.join(outpdf_dir, f"{new_file_name}.pdf")
    outpkl_path = os.path.join(outpkl_dir, f"{new_file_name}.pkl")
    outpkl_pT_corr_path = outpkl_path.replace(".pkl", "_pTcorrfactors.pkl")
    # See if these files already exist and overwrite if it's OK to do so.
    makeDirs(outpdf_dir, verbose)
    makeDirs(outpkl_dir, verbose)
    for f in [outpdf_path, outpkl_path, outpkl_pT_corr_path]:
        check_overwrite(f, overwrite)
    return (outpdf_path, outpkl_path, outpkl_pT_corr_path)

# def main():
assert not (do_binned_analysis and do_unbinned_analysis)
if not show_plots_to_screen: 
    r.gROOT.SetBatch(True)
# Prep your area.
outpdf_path, outpkl_path, outpkl_pT_corr_path = prep_files(eta_ls, pT_ls, outfile_prefix, outpdf_dir, outpkl_dir, 
                                                            regions, iters, overwrite=overwrite)
printer = CanvasPrinter(show_statsbox=False)
# Begin analysis.
muon_collection = MyMuonCollection()
muon_collection.extract_muons_from_Hmumu_file(infile_path, n_evts=max_n_evts, print_out_every=print_out_every,
                                    do_mu_pT_corr=False, 
                                    pT_corr_factor_dict=None,
                                    verbose=False)
# Sort muons into specified (eta, pT, qd0) bins.
muon_collection.sort_muons(eta_ls, pT_ls, verbose=verbose)
muon_collection.make_KinBin3Ds(regions=regions, 
                                min_muons_per_qd0_bin=min_muons_per_qd0_bin, 
                                verbose=verbose)
# Decide to do iter Gaus binning of qd0 or to do scatterplot of qd0.
if do_binned_analysis:
    muon_collection.do_3D_iter_gaus_fits(bins_dpTOverpT=bins_dpTOverpT, bins_qd0=bins_qd0, 
                                            x_lim_dpTOverpT=x_lim_dpTOverpT, x_lim_qd0=x_lim_qd0, 
                                            fit_whole_range_first_iter=fit_whole_range_first_iter,
                                            iters=iters,
                                            num_sigmas=num_sigmas,
                                            switch_to_binned_fit=switch_to_binned_fit,
                                            alarm_level=alarm_level,
                                            verbose=verbose)
    if make_dpTOverpT_vs_qd0_graphs:
        muon_collection.make_KinBin2D_graphs()
        muon_collection.collect_KinBin2D_plots()
        # Make a TMultiGraph for each eta bin.
        for eta_min,eta_max in zip(eta_ls[:-1], eta_ls[1:]):
            muon_collection.make_multigraph(eta_min, eta_max)
    muon_collection.collect_KinBin3D_plots()

    print(f"Drawing plots to:\n{outpdf_path}")
    canv = printer.canv
    canv.Print(outpdf_path + "[")
    # for eta_min,eta_max in zip(eta_ls[:-1], eta_ls[1:]):
    for mg in muon_collection.multigraph_ls:
        muon_collection.draw_mg_and_fits(mg)
        canv.Print(outpdf_path)
        r.gPad.Update()
    printer.make_plots_pretty(show_statsbox=True)
    # for gr in muon_collection.kinbin2d_plot_ls:
    #     gr.GetL
    #     color = kb2d.gr_dpTOverpT_vs_qd0.GetLineColor()
    #         kb2d.fit_line.SetLineColor()
    #         self.kinbin2d_plot_ls.extend(kb2d.gr_dpTOverpT_vs_qd0)

    # Now draw each multigraph's line, individually.
    for mg in muon_collection.multigraph_ls:
        gr_ls = list(mg.GetListOfGraphs())
        for gr in gr_ls:
            # color = ct + 1 if ct >= 5 else ct
            # gr.fit_line.SetLineColor(color)
            gr.Draw("ap")
            gr.fit_line.Draw("same")
            canv.Print(outpdf_path)
    # printer.draw_plots(muon_collection.kinbin2d_plot_ls, add_plots_to_pdf=True)
    printer.draw_plots(
        # muon_collection.kinbin2d_plot_ls + muon_collection.kinbin3d_plot_ls, 
        muon_collection.kinbin3d_plot_ls, 
        outpdf_path,
        add_plots_to_pdf=True)
    canv.Print(outpdf_path + "]")

elif do_unbinned_analysis:
    # Make a scatterplot of dpT/pT vs. qd0 values for each muon. 
    # Make all plots.
    muon_collection.make_inclusive_kinematic_plots()
    muon_collection.fill_KinBin_hists()
    muon_collection.make_KinBin2D_graphs()
    muon_collection.collect_all_plots()
    # Print plots into single PDF.
    print(f"Drawing all plots to:\n{outpdf_path}")
    printer.make_pdf_of_plots(muon_collection.plot_ls, outpdf_path)
    # muon_collection.save_to_pkl(muon_collection.pT_corr_factor_dict, outpkl_pT_corr_path)

if save_KinBin2D_pkl:
    muon_collection.overwrite_longlist_muon_info()
    muon_collection.save_to_pkl(muon_collection.KinBin2D_dict, outpkl_path)
if save_pT_corr_pkl:
    muon_collection.make_pT_corr_dict()
    muon_collection.save_to_pkl(muon_collection.pT_corr_factor_dict, outpkl_pT_corr_path)

# if __name__ == "__main__":
#     main()