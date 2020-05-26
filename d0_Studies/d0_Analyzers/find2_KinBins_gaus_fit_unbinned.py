"""
# PURPOSE: 
#   Load the dpT/pT data from data_dict.
#   Perform either a binned or unbinned fit over all those data. 
#   Plot the fits in a PDF. If unbinned, then make binned just for plotting. 
#   and save the stats in a readable '.stat' file and in a pickled dict.
# NOTES:   
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
# SYNTAX:  python script.py
# AUTHOR:  Jake Rosenzweig
# UPDATED: 2020-05-25
"""
import os
import sys
import math
import time
import pickle
import ROOT

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Local imports.
from d0_Studies.kinematic_bins import (equal_entry_bin_edges_eta_mod1)
from d0_Utils.d0_dicts import label_LaTeX_dict
from d0_Utils.d0_fns import make_binning_array, print_header_message
from d0_Utils.d0_cls import KinBin3D

from PyUtils.Utils_Physics import perc_diff
from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from PyUtils.Utils_Plotting import hist_y_label, make_1D_dist, ncolsrows_from_nplots, get_stats_1Dhist
from PyUtils.Utils_StatsAndFits import iterative_fit_gaus, iterative_fit_gaus_unbinned

#---------------------------#
#----- User Parameters -----#
#---------------------------#
# Load muon data.
inpath_data_dict_pkl = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/20200525_test02_smallskim_6regwith2000perreg__0p0_eta_2p3__100p0_pT_1000p0_GeV_datadict.pkl"
# inpath_data_dict_pkl = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/20200525_fullstats_6regwith2000perreg__0p0_eta_2p4__5p0_pT_1000p0_GeV_datadict.pkl"

# Outgoing dirs. 
outdir_plots = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Plots/hists_dpToverpT/"
outdir_kinbin_pkl = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/kinbins/"

overwrite = True
verbose = True
skip_bad_fit = True
do_unbinned_fit = True
iter_gaus = (True, 5)

# Kinematic to be plotted on all histograms. 
# Should not contain 1 or 2. 
# Acceptable values found in prepare_vaex_df().
kinem = "delta_pToverGenpT"  
x_bin_info = [-0.3, 0.3, 0.003]
x_zoom_range = [-0.35, 0.35]

# Cuts to make.
dR_max = 0.008
massZ_minmax_DY = [60, 120]
massZ_minmax_Jpsi = [2.9, 3.3]

# Leave this list empty if you want automatic detection.
# Will grab all sample_mods in data_dict.
# These are dict keys.
sample_mod_ls = [
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

# suffix   = "_{:.1f}_eta_{:.1f}".format(min(eta_ls), max(eta_ls))
# suffix  += "_{:.1f}_pT_{:.1f}_GeV".format(min(pT_ls), max(pT_ls))
# suffix += "_unbinned"
# suffix = make_str_title_friendly(suffix)

filename = os.path.split(inpath_data_dict_pkl)[1]
filename = filename.strip("_datadict.pkl")
if (do_unbinned_fit): 
    filename += "_unbinned"
stats_fullpath = os.path.join(outdir_plots, filename + "_kinbins.stat")
fullpath_kinbin_ls_pkl = os.path.join(outdir_kinbin_pkl, filename + "_kinbin_ls.pkl") 
fullpath_kinbin_dict_pkl = os.path.join(outdir_kinbin_pkl, filename + "_kinbin_dict.pkl") 

check_overwrite(stats_fullpath, overwrite=overwrite) 
check_overwrite(fullpath_kinbin_ls_pkl, overwrite=overwrite) 
check_overwrite(fullpath_kinbin_dict_pkl, overwrite=overwrite) 
# If file is there and you wish to overwrite, 
# make sure not to append to existing file; just get rid of it. 
if (overwrite):
    for file_ in [stats_fullpath, fullpath_kinbin_ls_pkl, fullpath_kinbin_dict_pkl]:
        try:
            os.remove(file_)
        except:
            pass

with open(stats_fullpath, "a") as f:
    f.write("#NOTE: hist_stats=[n_entries, mean, mean_err, std, std_err]\n\n")
    f.write("kinematic_variable, {}\n\n".format(kinem))

with open(inpath_data_dict_pkl, "rb") as f:
    data_dict = pickle.load(f)

# Unpack data.
eta_ls = data_dict["all_eta_bins"]
pT_ls = data_dict["all_pT_bins"]
if len(sample_mod_ls) == 0:
    sample_mod_ls = data_dict["sample_mod_ls"]

massZ_min_DY = massZ_minmax_DY[0]
massZ_max_DY = massZ_minmax_DY[1]
massZ_min_Jpsi = massZ_minmax_Jpsi[0]
massZ_max_Jpsi = massZ_minmax_Jpsi[1]

# plt.style.use('grid_multiple_plots')
ROOT.RooMsgService.instance().setStreamStatus(1,False)
plt.style.use("/home/rosedj1/.config/matplotlib/mpl_configdir/stylelib/grid_multiple_plots_python2.mplstyle")

x_bins, bin_width = make_binning_array(x_bin_info)

total_entries = 0
kinbin3D_ls_OLD = []  # Probably deprecate this.
kinbin_dict = {
    "all_eta_bins" : eta_ls,
    "all_pT_bins" : pT_ls,
    "sample_mod_ls" : sample_mod_ls,
    # "binninginfo_xbins_binwidth" : (x_bins, bin_width),
}

def get_grid_info(qd0_ls):
    """ Determine grid size per page. Return grid info."""
    n_plots_per_page = len(qd0_ls) - 1
    rows, cols = ncolsrows_from_nplots(n_plots_per_page)
    print("  This q*d0 bin edge list: {}".format(np.round(qd0_ls, decimals=5)))
    print("  This page will contain {} plots ({} x {} grid).\n".format(n_plots_per_page, rows, cols))
    return rows, cols

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

#----------------#
#----- Main -----#
#----------------#
for sample_mod in sample_mod_ls:
    # sample_mod is a key in data_dict.

    # Loop over eta regions.
    for k in range(len(eta_ls)-1):
        eta_min = eta_ls[k]
        eta_max = eta_ls[k+1]
        eta_range = [eta_min, eta_max]
        eta_key = f"eta_bin_left_edge={eta_min}"
        print(f"eta loop k={k}")

        kinbin_dict[eta_key] = {}
        # For each eta range, make a pdf. 
        filename_base = filename.split("__")[0]
        extra  = (
            f"__{eta_min:.1f}_eta_{eta_max:.1f}"
            f"__{min(pT_ls):.1f}_pT_{max(pT_ls):.1f}_GeV"
        )
        extra = make_str_title_friendly(extra)
        extra += ".pdf"

        outpath_pdf = os.path.join(outdir_plots, filename_base + extra)
        check_overwrite(outpath_pdf, overwrite=overwrite)    

        print(f"Making PDF {k+1}/{n_pdfs} for {sample_mod} sample...")
        t_start_page = time.perf_counter()
        with PdfPages(outpath_pdf) as pdf:
            # Within 1 eta region, scan the pT regions. 
            for j in range(len(pT_ls)-1):
                pT_min = pT_ls[j]
                pT_max = pT_ls[j+1]
                pT_range = [pT_min, pT_max]
                pT_key = f"pT_bin_left_edge={pT_min}"
                print("pT loop j,",j)

                kinbin_dict[eta_key][pT_key] = {}
                # Within 1 pT region, scan the q*d0 regions and make 1 page of PDF. 
                t_start = time.perf_counter()
                f = plt.figure()

                # Get the specific q*d0 list for this particular eta, pT region. 
                qd0_ls = [float(key.split("=")[1]) for key in data_dict[eta_key][pT_key][sample_mod].keys()]
                
                print("eta_key,",eta_key)
                print("pT_key,",pT_key)
                print("qd0_ls,",qd0_ls)

                rows, cols = get_grid_info(qd0_ls)

                # make a kinbin list for each (eta, pT) region.
                kinbin3D_ls = []

                for count in range(len(qd0_ls)-1):
                    ax = plt.subplot(rows,cols,count+1)
                    qd0_min = qd0_ls[count]
                    qd0_max = qd0_ls[count+1]
                    qd0_range = [qd0_min, qd0_max]
                    qd0_key = f"qd0_bin_left_edge={qd0_min}"
                    print(f"qd0 loop: {count}")
                    print(f"qd0_range={qd0_range}")
                    
                    x_label = label_LaTeX_dict[kinem]["independent_label"]
                    x_units = label_LaTeX_dict[kinem]["units"]
                    y_label = hist_y_label(bin_width, x_units)
                    if len(x_units) > 0:
                        x_label += " [{}]".format(x_units)

                    cuts  = r"$%.2f < \left| \eta^{\mathrm{REC}} \right| < %.2f$," % (eta_min, eta_max) + "\n"
                    cuts += r"$%.1f <$ %s $< %.1f$ GeV," % (pT_min, p_str_latex, pT_max) + "\n"
                    cuts += r"$%.4f < q(\mu)*d_{0}^{\mathrm{%s}} < %.4f$, " % (qd0_min, "BS", qd0_max) + "\n"
                    cuts += r"$\Delta R < %.3f$, " % (dR_max)
                    if "DY" in sample_mod:
                        cuts += "\n" + r"$%.1f < m_{Z} < %.1f$ GeV,  " % (massZ_min_DY, massZ_max_DY)
                    if "Jpsi" in sample_mod:
                        cuts += "\n" + r"$%.1f < m_{J/\psi} < %.1f$ GeV" % (massZ_min_Jpsi, massZ_max_Jpsi)
                    print("cuts str:\n",cuts)
                    
                    dpToverpT_vals = data_dict[eta_key][pT_key][sample_mod][qd0_key]['dpToverpT_vals']
                    stats_ls_pT    = data_dict[eta_key][pT_key][sample_mod][qd0_key]['stats_ls_pT']
                    stats_ls_qd0   = data_dict[eta_key][pT_key][sample_mod][qd0_key]['stats_ls_qd0']

                    # Make sure that everyone agrees on the number of muons:
                    assert stats_ls_pT[0] == stats_ls_qd0[0]
                    assert stats_ls_pT[0] == len(dpToverpT_vals)

                    if (do_unbinned_fit): print_header_message("--- Doing unbinned fit ---")
                    # Plot the kinem hist.
                    ax, bin_vals, bin_edges, stats_binned = make_1D_dist(
                                                        ax=ax, 
                                                        data=dpToverpT_vals,
                                                        x_limits=x_zoom_range,
                                                        x_bins=x_bins, 
                                                        x_label=x_label, 
                                                        y_label=y_label,
                                                        title="",
                                                        y_max=-1,
                                                        log_scale=False)
                    print(f"stats_binned ls: {stats_binned}")
                    ax.text(0.025, 0.78, cuts, horizontalalignment='left', verticalalignment='center', 
                    transform=ax.transAxes, bbox=dict(boxstyle='square', facecolor='white', alpha=0.9))

                    if (iter_gaus[0]):
                        # Do iterative fitting procedure.
                        num_iters = iter_gaus[1]
                        # Produce a dictionary of stats for these fits.

                        if (do_unbinned_fit):
                            # Returns a dictionary with different keys from the binned fit. 
                            fit_stats_dict, ax = iterative_fit_gaus_unbinned(num_iters, dpToverpT_vals,
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

                    with open(stats_fullpath, "a") as f:
                        f.write(f"range_eta,  {eta_range}\n")
                        f.write(f"range_pT,   {pT_range}\n")
                        f.write(f"range_qd0,  {qd0_range}\n")
                        f.write(f"hist_stats, {stats_binned}\n")
                        f.write(f"fit_stats:\n")
                        for key,val in fit_stats_dict.items():
                            f.write(f"{key}, {val}\n")
                        f.write("\n")
                        f.close()

                # End qd0 loop.
                kinbin_dict[eta_key][pT_key][sample_mod] = kinbin3D_ls

                t_end = time.perf_counter()
                plt.tight_layout()
                pdf.savefig()
                plt.close("all")
                
                msg = f"Page {j+1}/{n_pages} made for {sample_mod} sample. Time taken: {(t_end - t_start):.2f} s"
                print_header_message(msg)
                print()

            # End pT loop. Make next page.
        # PDF made.
        t_end_page = time.perf_counter()
        msg_made_pdf = f"PDF {k+1}/{n_pdfs} made"
        print_header_message(msg_made_pdf, pad_char="@", n_center_pad_chars=5)
        print(f"location:\n  {outpath_pdf}")
        print(f"(took {(t_end_page - t_start_page):.2f} s)\n")
        # Save this 1 eta reg, 1 pT reg, and all q*d0 plots on one page. 

    # Go to next pT reg and next page.
    print("All PDFs created.")
# End sample loop.

with open(fullpath_kinbin_ls_pkl,'wb') as output:
    pickle.dump(kinbin3D_ls, output, protocol=2)
print("[INFO] KinBin3D list written to pickle file:\n{}\n".format(fullpath_kinbin_ls_pkl))

with open(fullpath_kinbin_dict_pkl,'wb') as output:
    pickle.dump(kinbin_dict, output, protocol=2)
print("[INFO] KinBin3D dict written to pickle file:\n{}\n".format(fullpath_kinbin_dict_pkl))
    
# total_muons_original = vdf_concat_MC_2017_DY.count() + vdf_concat_MC_2017_Jpsi.count()
# print("Total muons: {}".format(total_muons_original))
# perdif = perc_diff(total_entries, total_muons_original)
# print("Total muons found across all bins: {} (perc. diff. = {:.2f}%)".format(total_entries, perdif))