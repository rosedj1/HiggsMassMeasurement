"""
# PURPOSE: 
#   Save the delta_pT/pT data for each (eta, pT, q*d0) cube. 
#   Also saves the statistics of pT and q*d0 (like the avg, std, and errors)
#   into a dict, which is stored in a '.pkl'.
#   All these data are to be eventually fed into a RooFit unbinned fit. 
#   Load the dpT/pT data from data_dict.
#   Perform either a binned or unbinned fit over all those data. 
#   Plot the fits in a PDF. If unbinned, then make binned just for plotting. 
#   and save the stats in a readable '.stat' file and in a pickled dict.
# NOTES:   
#   User must specify the input pkl file which contains the equal-entry bin edges:
#     pickled_dict = {eta_bin_left_edge : {
                        pT_bin_left_edge : {
                            "equalentry_qd0ls_DY" : [],
                            "equalentry_qd0ls_Jpsi" : [],
                            "equalentry_qd0ls_DY+Jpsi" : [],
                        }
                    }
#   THANKS TO CONDA VIRTUAL ENVIRONMENT, WE CAN RUN THIS IN PYTHON 3!!!
#   Performs unbinned Gaussian fits of distributions, can show the fit stats
#   on the plots, and saves the fit info in a ".stat" file.
#   Optional: make a pickle file of the list of KinBin3D objects. 
#   One PDF is made per eta bin. 
#   Each page of the PDF represents a different pT bin.
#   Many plots are shown on a single page; these are different q*d0 bins.
#   Each page can have a various number of q*d0 In total, N q*d0BS distributions are made per page.
#   User must specify the dict which contains the bins:
#              pickled_dict = {eta_bin_min : {pT_bin_min : [qd0_ls]} }
#   N is determined automatically based on each len(qd0_ls).
#   User should check User Parameters.
#   On a full run of 13 eta bins, 12 pT bins, and 20 q*d0 regions (so 3120 "cubes"), 
#     time elapsed ~6 hrs. This is ~8.5 plots/min.
#   Binning is automatically detected.
#   User should check User Parameters.
# SYNTAX:  python script.py
# AUTHOR:  Jake Rosenzweig
# UPDATED: 2020-06-24
"""
import time
t_prog_start = time.perf_counter()
import os
import sys
import math
import pickle
import ROOT
# import argparse 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Local imports.
from Utils_vaex.vaex_fns import prepare_vaex_df, vaex_apply_masks
from Samples.sample_info import Sample
from d0_Utils.d0_dicts import label_LaTeX_dict
from d0_Utils.d0_fns import make_binning_array, print_header_message
from d0_Utils.d0_cls import KinBin3D

from Utils_Python.Utils_Physics import perc_diff
from Utils_Python.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from Utils_Python.Utils_Plotting import hist_y_label, make_1D_dist, ncolsrows_from_nplots, get_stats_1Dhist
from Utils_Python.Utils_StatsAndFits import iterative_fit_gaus, iterative_fit_gaus_unbinned
# def ParseOption():
#     parser = argparse.ArgumentParser(description='submit all')
#     parser.add_argument('--pklfilename', dest='filename_base', type=str, help='') 
#     parser.add_argument('--verbose', dest='verbose', type=int, default=1, help='')
#     parser.add_argument('--overwrite', dest='overwrite', type=int, default=0, help='')  
    
#     args = parser.parse_args()                                                                                         
#     return args          
                                                                                                         
# args = ParseOption()          

#---------------------------#
#----- User Parameters -----#
#---------------------------#
# Dictionary which contains equal-entry q*d0 bin edges.
# inpath_equalentry_pickl_dict = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/20200526_fullstats_BA_OL_EC_12regwith3000perreg__0p0_eta_2p4__5p0_pT_200p0_GeV.pkl"
# inpath_equalentry_pickl_dict = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/20200531_fullstats_sameasFilippodeltaRcut_12regwith3000perreg__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"
# inpath_equalentry_pickl_dict = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/20200624_MC2018JpsiDY_test01_6regwith3000perreg__0p0_eta_0p4__5p0_pT_1000p0_GeV.pkl"
inpath_equalentry_pickl_dict = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/20200624_fullstats_MC2018JpsiDY_12regwith3000perreg__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"

# Outgoing dirs. 
outdir_plots = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Plots/hists_dpToverpT/"
outdir_kinbin_pkl = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/kinbins/"

MC_2018_Jpsi_hdf5 = Sample("MC", "2018", "Jpsi", "hdf5")
MC_2018_DY_hdf5 = Sample("MC", "2018", "DY", "hdf5")

# Filename of outdict is automatic.
overwrite = False
verbose = True
skip_bad_fit = True
do_unbinned_fit = True
iter_gaus = (True, 5)

# Kinematic to be plotted on all histograms. 
# Should not contain 1 or 2. 
# Acceptable values found in prepare_vaex_df().
kinem = "delta_pToverGenpT"  
x_bin_info = [-0.35, 0.35, 0.0025]
x_zoom_range = [-0.35, 0.35]

# Selections per sample.
# dR_max = 0.008
dR_max_DY = 0.002
dR_max_Jpsi = 0.005
massZ_minmax_DY = [60, 120]
massZ_minmax_Jpsi = [2.9, 3.3]

# Choose the samples you want to run over.
# Leave this list empty if you want to run over 
# all available samples found in equalentry_pickl_dict
# These are dict keys.
sample_name_ls = [
    # "DY",
    # "Jpsi",
    "DY+Jpsi"
]

p_str_latex = "$p_{T}$"
#----------------------#
#----- Automatons -----#
#----------------------#
# Prep directories and file names. 
makeDirs(outdir_plots)
makeDirs(outdir_kinbin_pkl)

filename = os.path.split(inpath_equalentry_pickl_dict)[1]
filename = filename.strip(".pkl")
if (do_unbinned_fit): 
    filename += "_unbinned"
fullpath_stats = os.path.join(outdir_plots, filename + "_kinbins.stat")
fullpath_kinbin_ls_pkl = os.path.join(outdir_kinbin_pkl, filename + "_kinbin_ls.pkl") 
fullpath_kinbin_dict_pkl = os.path.join(outdir_kinbin_pkl, filename + "_kinbin_dict.pkl") 

check_overwrite(fullpath_stats, overwrite=overwrite) 
check_overwrite(fullpath_kinbin_ls_pkl, overwrite=overwrite) 
check_overwrite(fullpath_kinbin_dict_pkl, overwrite=overwrite) 
# If file is there and you wish to overwrite, 
# make sure not to append to existing file; just get rid of it. 
if (overwrite):
    for file_ in [fullpath_stats, fullpath_kinbin_ls_pkl, fullpath_kinbin_dict_pkl]:
        try:
            os.remove(file_)
        except:
            pass

with open(fullpath_stats, "a") as f:
    f.write("#NOTE: hist_stats=[n_entries, mean, mean_err, std, std_err]\n\n")
    f.write("kinematic_variable, {}\n\n".format(kinem))

with open(inpath_equalentry_pickl_dict, "rb") as f:
    equalentry_binedge_dict = pickle.load(f)

# Unpack data.
# Binning automatically detected from input pickle file.
eta_ls = equalentry_binedge_dict["all_eta_bins"]
pT_ls = equalentry_binedge_dict["all_pT_bins"]
if len(sample_name_ls) == 0:
    # Automatically detect all samples in dict.
    sample_name_ls = equalentry_binedge_dict["sample_name_ls"]

massZ_min_DY = massZ_minmax_DY[0]
massZ_max_DY = massZ_minmax_DY[1]
massZ_min_Jpsi = massZ_minmax_Jpsi[0]
massZ_max_Jpsi = massZ_minmax_Jpsi[1]

total_num_muons_DY = 0
total_num_muons_Jpsi = 0
total_entries = 0
kinbin3D_ls_OLD = []  # Probably deprecate this.
kinbin_dict = {
    "all_eta_bins" : eta_ls,
    "all_pT_bins" : pT_ls,
    "sample_name_ls" : sample_name_ls,
    # "binninginfo_xbins_binwidth" : (x_bin_edges, bin_width),
}

#---------------------#
#--- Plot Settings ---#
#---------------------#
ROOT.RooMsgService.instance().setStreamStatus(1,False)

x_bin_edges, bin_width = make_binning_array(x_bin_info)

def get_grid_info(qd0_ls):
    """ Determine grid size per page. Return grid info."""
    n_plots_per_page = len(qd0_ls) - 1
    rows, cols = ncolsrows_from_nplots(n_plots_per_page)
    print("  This q*d0 bin edge list: {}".format(np.round(qd0_ls, decimals=5)))
    print("  This page will contain {} plots ({} x {} grid).\n".format(n_plots_per_page, rows, cols))
    return rows, cols

x_label = label_LaTeX_dict[kinem]["independent_label"]
x_units = label_LaTeX_dict[kinem]["units"]
y_label = hist_y_label(bin_width, x_units)
if len(x_units) > 0:
    x_label += " [{}]".format(x_units)

# Determine PDF and page info.
n_pdfs = len(eta_ls)-1
n_pages = len(pT_ls) - 1
msg = "Making {} PDF.".format(n_pdfs)
if n_pdfs > 1:
    msg = msg.replace("PDF", "PDFs")
print(msg)
print("Each PDF will contain {} pages.\n".format(n_pages))
print("|eta| regions: {}\n".format(np.round(eta_ls, decimals=2)))
print("pT regions: {}\n".format(np.round(pT_ls, decimals=2)))

print(
    f"eta_ls: {eta_ls}\n",
    f"pT_ls:  {pT_ls}",
)

#---------------------------------#
#----- Unpack kinematic data -----#
#---------------------------------#
eta_arr_DY     = MC_2018_DY_hdf5.vdf_prepped.evaluate("eta")
eta_arr_Jpsi   = MC_2018_Jpsi_hdf5.vdf_prepped.evaluate("eta")
pT_arr_DY      = MC_2018_DY_hdf5.vdf_prepped.evaluate("pT")
pT_arr_Jpsi    = MC_2018_Jpsi_hdf5.vdf_prepped.evaluate("pT")
qd0_arr_DY     = MC_2018_DY_hdf5.vdf_prepped.evaluate("qd0BS")
qd0_arr_Jpsi   = MC_2018_Jpsi_hdf5.vdf_prepped.evaluate("qd0BS")
massZ_arr_DY   = MC_2018_DY_hdf5.vdf_prepped.evaluate("massZ")
massZ_arr_Jpsi = MC_2018_Jpsi_hdf5.vdf_prepped.evaluate("massZ")
dR_arr_DY      = MC_2018_DY_hdf5.vdf_prepped.evaluate("delta_R")
dR_arr_Jpsi    = MC_2018_Jpsi_hdf5.vdf_prepped.evaluate("delta_R")

# Prep the masks that don't change. 
mask_massZ_DY   = (massZ_min_DY < massZ_arr_DY) & (massZ_arr_DY < massZ_max_DY)
mask_massZ_Jpsi = (massZ_min_Jpsi < massZ_arr_Jpsi) & (massZ_arr_Jpsi < massZ_max_Jpsi)
mask_dR_DY   = (dR_arr_DY < dR_max_DY)
mask_dR_Jpsi = (dR_arr_Jpsi < dR_max_Jpsi)

# Kinematic to be plotted in each KinBin3D.
# Should not contain 1 or 2. 
# Acceptable values found in prepare_vaex_df().
kinem_arr_DY_dpToverGenpT = MC_2018_DY_hdf5.vdf_prepped.evaluate("delta_pToverGenpT")
kinem_arr_DY_pT = MC_2018_DY_hdf5.vdf_prepped.evaluate("pT")
kinem_arr_DY_qd0 =  MC_2018_DY_hdf5.vdf_prepped.evaluate("qd0BS")
kinem_arr_Jpsi_dpToverGenpT = MC_2018_Jpsi_hdf5.vdf_prepped.evaluate("delta_pToverGenpT")
kinem_arr_Jpsi_pT = MC_2018_Jpsi_hdf5.vdf_prepped.evaluate("pT")
kinem_arr_Jpsi_qd0 =  MC_2018_Jpsi_hdf5.vdf_prepped.evaluate("qd0BS")

#----------------#
#----- Main -----#
#----------------#
for sample_name in sample_name_ls:
    # sample_name is a key in kinbin_dict.
    # Loop over eta regions.
    for k in range(len(eta_ls)-1):
        eta_min = eta_ls[k]
        eta_max = eta_ls[k+1]
        eta_range = [eta_min, eta_max]
        eta_key = f"eta_bin={eta_range}"
        print(f"eta loop k={k}")

        kinbin_dict[eta_key] = {}
        # For each eta range, make a pdf. 
        filename_base = filename.split("__")[0]
        extra  = (
            f"_{sample_name}2018"
            f"_{eta_min}_eta_{eta_max}"
            f"_{min(pT_ls):.1f}_pT_{max(pT_ls):.1f}_GeV"
        )
        extra = make_str_title_friendly(extra)
        extra += ".pdf"

        outpath_pdf = os.path.join(outdir_plots, filename_base + extra)
        check_overwrite(outpath_pdf, overwrite=overwrite)    

        print(f"Making PDF {k+1}/{n_pdfs} for {sample_name} sample...")
        t_start_page = time.perf_counter()
        with PdfPages(outpath_pdf) as pdf:
            # Within 1 eta region, scan the pT regions. 
            for j in range(len(pT_ls)-1):
                pT_min = pT_ls[j]
                pT_max = pT_ls[j+1]
                pT_range = [pT_min, pT_max]
                pT_key = f"pT_bin={pT_range}"
                print(f"pT loop j={j}")

                # Within 1 pT region, scan the q*d0 regions and make 1 page of PDF. 
                t_start = time.perf_counter()
                f = plt.figure()

                # Get the equal-entry qd0 bin edges for this sample_name.
                qd0_ls_key = f"equalentry_qd0ls_{sample_name}"
                qd0_ls = equalentry_binedge_dict[eta_key][pT_key][qd0_ls_key]

                # Save this qd0_ls.
                kinbin_dict[eta_key][pT_key] = {}
                kinbin_dict[eta_key][pT_key][sample_name] = {}
                # kinbin_dict[eta_key][pT_key][sample_name] is a dictionary in which
                # this key is unlike the others. 
                kinbin_dict[eta_key][pT_key][sample_name][qd0_ls_key] = qd0_ls
                print("eta_key,",eta_key)
                print("pT_key,",pT_key)
                print("qd0_ls,",qd0_ls)

                rows, cols = get_grid_info(qd0_ls)

                if ((rows * cols) > 6) and (iter_gaus[1] > 3):
                    # Lots of fit info will be printed to legend. Shrink font. 
                    plt.style.use("/home/rosedj1/.config/matplotlib/mpl_configdir/stylelib/grid_morethan6plots_manyiters.mplstyle")
                else:
                    plt.style.use("/home/rosedj1/.config/matplotlib/mpl_configdir/stylelib/grid_multiple_plots.mplstyle")

                # Make a kinbin list for this (eta, pT) region.
                kinbin3D_ls = []
                for count in range(len(qd0_ls)-1):
                    # Go through each qd0_range and make 1 KinBin3D obj.
                    ax = plt.subplot(rows,cols,count+1)
                    qd0_min = qd0_ls[count]
                    qd0_max = qd0_ls[count+1]
                    qd0_range = [qd0_min, qd0_max]
                    qd0_key = f"qd0_bin={qd0_range}"
                    print(f"qd0 loop: {count}")
                    print(f"qd0_range={qd0_range}")
                    
                    kinbin_dict[eta_key][pT_key][sample_name][qd0_key] = {}

                    if (verbose):
                        print(f"\nqd0_loop={count}")
                        print(f"sample_name={sample_name}")
                        print(f"eta_range={eta_range}")
                        print(f"pT_range={pT_range}")
                        print(f"qd0_range={qd0_range}")

                    # Make masks, apply cuts, get data.
                    try:
                        # Prepare masks.
                        # These must be labelled with "DY" since they are of DY's size.
                        if "DY" in sample_name:
                            print(f"Applying DY masks")
                            # Apply masks to DY sample.
                            mask_eta_DY = (eta_min < np.abs(eta_arr_DY)) & (np.abs(eta_arr_DY) < eta_max)
                            mask_pT_DY = (pT_min < pT_arr_DY) & (pT_arr_DY < pT_max)
                            mask_qd0_DY = (qd0_min < qd0_arr_DY) & (qd0_arr_DY < qd0_max)
                            # Combine masks.
                            all_masks_DY = mask_eta_DY & mask_pT_DY & mask_qd0_DY & mask_massZ_DY & mask_dR_DY
                            # Apply masks.
                            data_DY_dpToverGenpT = kinem_arr_DY_dpToverGenpT[all_masks_DY]
                            data_DY_pT = kinem_arr_DY_pT[all_masks_DY]
                            data_DY_qd0 = kinem_arr_DY_qd0[all_masks_DY]

                        if "Jpsi" in sample_name:
                            print(f"Applying Jpsi masks")
                            # Apply masks to J/psi sample.
                            mask_eta_Jpsi = (eta_min < np.abs(eta_arr_Jpsi)) & (np.abs(eta_arr_Jpsi) < eta_max)
                            mask_pT_Jpsi = (pT_min < pT_arr_Jpsi) & (pT_arr_Jpsi < pT_max)
                            mask_qd0_Jpsi = (qd0_min < qd0_arr_Jpsi) & (qd0_arr_Jpsi < qd0_max)
                            # Combine masks.
                            all_masks_Jpsi = mask_eta_Jpsi & mask_pT_Jpsi & mask_qd0_Jpsi & mask_massZ_Jpsi & mask_dR_Jpsi
                            # Apply masks.
                            data_Jpsi_dpToverGenpT = kinem_arr_Jpsi_dpToverGenpT[all_masks_Jpsi]
                            data_Jpsi_pT = kinem_arr_Jpsi_pT[all_masks_Jpsi]
                            data_Jpsi_qd0 = kinem_arr_Jpsi_qd0[all_masks_Jpsi]

                    except TypeError:
                        print(f"[WARNING] sample_name={sample_name} threw a TypeError.")
                        print("There were probably 0 muons found in this region. Continuing anyway.")
                        # kinbin_dict[eta_key][pT_key][sample_name][qd0_key]["dpToverpT_vals"] = [None]
                        kinbin_dict[eta_key][pT_key][sample_name][qd0_key]["stats_ls_pT"] = [None]
                        kinbin_dict[eta_key][pT_key][sample_name][qd0_key]["stats_ls_qd0"] = [None]
                        continue

                    # Analyze the data with selections.
                    if sample_name == "DY":   
                        mask = all_masks_DY
                        data = data_DY_dpToverGenpT
                        stats_ls_pT = get_stats_1Dhist(data_DY_pT)
                        stats_ls_qd0 = get_stats_1Dhist(data_DY_qd0)

                        n_DY = len(data)
                        total_num_muons_DY += n_DY
                        if (verbose): 
                            print(f"  Number of DY muons in this cube: {n_DY}")
                            print(f"  Cumulative number of DY muons:   {total_num_muons_DY}")
                    elif sample_name == "Jpsi": 
                        mask = all_masks_Jpsi
                        data = data_Jpsi_dpToverGenpT
                        stats_ls_pT = get_stats_1Dhist(data_Jpsi_pT)
                        stats_ls_qd0 = get_stats_1Dhist(data_Jpsi_qd0)

                        n_Jpsi = len(data)
                        total_num_muons_Jpsi += n_Jpsi
                        if (verbose): 
                            print(f"  Number of Jpsi muons in this cube: {n_Jpsi}")
                            print(f"  Cumulative number of Jpsi muons:   {total_num_muons_Jpsi}")
                    elif sample_name == "DY+Jpsi":
                        print("Concatenating DY+Jpsi data together")
                        mask = np.concatenate((all_masks_DY, all_masks_Jpsi))
                        data = np.concatenate((data_DY_dpToverGenpT, data_Jpsi_dpToverGenpT))
                        stats_ls_pT = get_stats_1Dhist(np.concatenate((data_DY_pT, data_Jpsi_pT)))
                        stats_ls_qd0 = get_stats_1Dhist(np.concatenate((data_DY_qd0, data_Jpsi_qd0)))
                    else: 
                        raise KeyError(f"sample_name={sample_name} not recognized.")
                        
                    # A couple of checks.
                    num_passed = len(data)
                    sum_mask = np.sum(mask)
                    assert sum_mask == num_passed

                    # Data acquired. Fit and plot it.
                    cuts  = r"$%.2f < \left| \eta^{\mathrm{REC}} \right| < %.2f$," % (eta_min, eta_max) + "\n"
                    cuts += r"$%.1f <$ %s $< %.1f$ GeV," % (pT_min, p_str_latex, pT_max) + "\n"
                    cuts += r"$%.4f < q(\mu)*d_{0}^{\mathrm{%s}} < %.4f$, " % (qd0_min, "BS", qd0_max) + "\n"
                    if "Jpsi" in sample_name:
                        cuts += "\n" + r"$%.1f < m_{J/\psi} < %.1f$ GeV,  $\Delta R < %.3f$" % (massZ_min_Jpsi, massZ_max_Jpsi, dR_max_Jpsi)
                    if "DY" in sample_name:
                        cuts += "\n" + r"$%.1f < m_{Z} < %.1f$ GeV,  $\Delta R < %.3f$" % (massZ_min_DY, massZ_max_DY, dR_max_DY)
                    if (verbose): print("cuts str:\n",cuts)

                    # Make sure that everyone agrees on the number of muons:
                    assert stats_ls_pT[0] == stats_ls_qd0[0]
                    assert stats_ls_pT[0] == len(data)
    
                    if (do_unbinned_fit): print_header_message("--- Doing unbinned fit ---")
                    # Plot the kinem hist.
                    ax, bin_vals, bin_edges, stats_binned = make_1D_dist(
                                                        ax=ax, 
                                                        data=data,
                                                        x_limits=x_zoom_range,
                                                        x_bin_edges=x_bin_edges, 
                                                        x_label=x_label, 
                                                        y_label=y_label,
                                                        title="",
                                                        y_max=-1,
                                                        log_scale=False, color=None, leg_loc=None, display="sci")
                    print(f"stats_binned ls: {stats_binned}")
                    ax.text(0.025, 0.78, cuts, horizontalalignment='left', verticalalignment='center', 
                    transform=ax.transAxes, bbox=dict(boxstyle='square', facecolor='white', alpha=0.9))

                    # Optional iterative Gaussian fits.
                    if (iter_gaus[0]):
                        # Do iterative fitting procedure.
                        num_iters = iter_gaus[1]
                        # Produce a dictionary of stats for these fits.

                        if (do_unbinned_fit):
                            # Returns a dictionary with different keys from the binned fit. 
                            fit_stats_dict, ax = iterative_fit_gaus_unbinned(num_iters, data,
                                                                            bin_edges=bin_edges, 
                                                                            bin_vals=bin_vals,
                                                                            num_sigmas=2,
                                                                            ax=ax, draw_on_axes=True, verbose=verbose)
                        else: 
                            # Do binned fit.
                            # The Gaus coeff is usually around the max height of Gaus curve.
                            best_guess_coeff = np.max(bin_vals)  
                            best_guess_mean = stats_binned[1]
                            best_guess_stdev = stats_binned[3]

                            fit_stats_dict, ax = iterative_fit_gaus(num_iters, 
                                                                bin_edges, bin_vals, 
                                                            param_guess=[best_guess_coeff, 
                                                                        best_guess_mean, 
                                                                        best_guess_stdev],
                                                            param_bounds=([0,-1000,-100], [999999999,1000,100]),
                                                            ax=ax, draw_on_axes=True, verbose=False, skip_bad_fit=skip_bad_fit)
                    else:
                        fit_stats_dict = None
                    
                    # Save all them tasty fit statistics. 
                    # Note that the fit_stats_dict will have different keys,
                    # depending on which fit was done (binned vs. unbinned).
                    fit_type = "unbinned" if do_unbinned_fit else "binned"
                    cube = KinBin3D(eta_range=eta_range,
                                    pT_range=pT_range, 
                                    qd0_range=qd0_range,
                                    n_entries=stats_ls_pT[0],
                                    kinem=kinem,
                                    fit_stats_dict=fit_stats_dict,
                                    fit_type=fit_type,
                                    pT_stats_ls=stats_ls_pT,
                                    qd0_stats_ls=stats_ls_qd0,
                                    cut_str=cuts
                            )
                    kinbin3D_ls_OLD.append(cube)  # Probably deprecate. Big, unorganized list.
                    kinbin3D_ls.append(cube)

                    kinbin_dict[eta_key][pT_key][sample_name][qd0_key]["stats_ls_pT"] = stats_ls_pT
                    kinbin_dict[eta_key][pT_key][sample_name][qd0_key]["stats_ls_qd0"] = stats_ls_qd0

                    with open(fullpath_stats, "a") as f:
                        f.write(f"range_eta,  {eta_range}\n")
                        f.write(f"range_pT,   {pT_range}\n")
                        f.write(f"range_qd0,  {qd0_range}\n")
                        f.write(f"hist_stats, {stats_binned}\n")
                        f.write(f"fit_stats:\n")
                        for key,val in fit_stats_dict.items():
                            f.write(f"{key}, {val}\n")
                        f.write("\n")
                        f.close()
                    # Finished this q*d0 bin. 
                # End qd0 loop.
                kinbin_dict[eta_key][pT_key][sample_name]["kinbin3D_ls"] = kinbin3D_ls

                t_end = time.perf_counter()
                plt.tight_layout()
                pdf.savefig()
                plt.close("all")
                
                msg = f"Page {j+1}/{n_pages} made for {sample_name} sample. Time taken: {(t_end - t_start):.2f} s"
                print_header_message(msg)
                print()
                # End all q*d0 regions for this list.
                # Save this 1 eta reg, 1 pT reg, and all q*d0 plots on one page. 
            # End pT loop. Make next page.
        # PDF made and closed.
        t_end_page = time.perf_counter()
        msg_made_pdf = f"PDF {k+1}/{n_pdfs} made"
        print_header_message(msg_made_pdf, pad_char="@", n_center_pad_chars=5)
        print(f"location:\n  {outpath_pdf}")
        print(f"(took {(t_end_page - t_start_page):.2f} s)\n")
    # End eta loop.
# End sample loop.
print("All PDFs created.")

with open(fullpath_kinbin_ls_pkl,'wb') as output:
    pickle.dump(kinbin3D_ls_OLD, output, protocol=2)
print("[INFO] KinBin3D list written to pickle file:\n{}\n".format(fullpath_kinbin_ls_pkl))

with open(fullpath_kinbin_dict_pkl,'wb') as output:
    pickle.dump(kinbin_dict, output, protocol=2)
print("[INFO] KinBin3D dict written to pickle file:\n{}\n".format(fullpath_kinbin_dict_pkl))

t_prog_end = time.perf_counter()
print(f"[INFO] Total time taken: {t_prog_end - t_prog_start} sec.")
# total_muons_original = MC_2018_DY_hdf5.vdf_prepped.count() + MC_2018_Jpsi_hdf5.vdf_prepped.count()
# total_muons_found = total_num_muons_DY + total_num_muons_Jpsi
# perdif = perc_diff(total_muons_found, total_muons_original)
# print(f"Total muons in files: {total_muons_original}")
# print(f"Total muons found across all bins: {total_muons_found} (perc. diff. = {perdif:.2f}%)")