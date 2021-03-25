# PURPOSE: Make PDFs of any kinematic distribution (e.g., delta_pT/pT), 
#          in specified eta, pT, and q*d0 bins.
#          Can perform Gaussian fits of distributions and put the fit stats
#          on the plots. 
# NOTES:   One PDF is made per eta bin. 
#          Each page of the PDF represents a different pT bin.
#          Many plots are shown on a single page; these are different q*d0 bins.
#          In total, N q*d0BS distributions are made per page.
#          N is determined automatically based on len(qd0_ls).
#          User should check User parameters.
# SYNTAX:  python <script>.py
# AUTHOR:  Jake Rosenzweig
# UPDATED: 2020-05-14

import os
import sys
import math
import time

import numpy as np
import matplotlib.pyplot as plt

# Local imports.
from vaex_Utils.vaex_dataframes import (vdf_MC_2017_DY, vdf_MC_2017_Jpsi, vdf_MC_2017_DY,
                                        prepare_vaex_df, vaex_apply_masks)
from d0_Studies.kinematic_bins import (equal_entry_bin_edges_eta_mod1)
from d0_Utils.d0_dicts import label_LaTeX_dict
from d0_Utils.d0_fns import make_binning_array

from Utils_Python.printing import print_header_message
from PyUtils.Utils_Physics import perc_diff
from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from PyUtils.Utils_Plotting import hist_y_label, make_1D_dist, ncolsrows_from_nplots
from PyUtils.Utils_StatsAndFits import iterative_fit_gaus

from matplotlib.backends.backend_pdf import PdfPages

#---------------------------#
#----- User Parameters -----#
#---------------------------#
# Samples:
vdf_concat_MC_2017_DY = prepare_vaex_df(vdf_MC_2017_DY)
vdf_concat_MC_2017_Jpsi = prepare_vaex_df(vdf_MC_2017_Jpsi)
outdir = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/hists_dpToverpT/MC/2017/"
filename_base = "Test08_overlapregion_MC2017DYandJpsi_wholenumber_pT__smallestRMS_qd0"
overwrite = False
verbose = True

# Binning.
eta_ls = equal_entry_bin_edges_eta_mod1_wholenum[2:4]
pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum[8:]
qd0_ls = binedges_qd0_equalentry_smallest_qd0RMS[4:9]

# Should not contain 1 or 2. 
# Acceptable values found in label_LaTeX_dict.
kinem = "delta_pToverGenpT"  
x_bin_info = [-0.1, 0.1, 0.002]
x_zoom_range = [-0.15, 0.15]
iter_gaus = (True, 3)

# Cuts to make.
dR_max = 0.008
massZ_minmax_DY = [60, 120]
massZ_minmax_Jpsi = [2.9, 3.3]

p_str_latex = "$p_{T}$"
#----------------------#
#----- Automatons -----#
#----------------------#
massZ_min_DY = massZ_minmax_DY[0]
massZ_max_DY = massZ_minmax_DY[1]
massZ_min_Jpsi = massZ_minmax_Jpsi[0]
massZ_max_Jpsi = massZ_minmax_Jpsi[1]

# Prep directories and file names. 
makeDirs(outdir)

stats_filename = filename_base
extra_   = "__{}_eta_{}".format(min(eta_ls), max(eta_ls))
extra_  += "__{:.1f}_pT_{:.1f}_GeV".format(min(pT_ls), max(pT_ls))
extra_ = make_str_title_friendly(extra_)
extra_ += ".stat"

stats_fullpath = os.path.join(outdir, stats_filename + extra_)
check_overwrite(stats_fullpath, overwrite=overwrite) 
# If file is there and you wish to overwrite, 
# make sure not to append to existing file; just get rid of it. 
if (overwrite):
    try:
        os.remove(stats_fullpath)
    except:
        pass

with open(stats_fullpath, "a") as f:
    f.write("#NOTE: hist_stats=[n_entries, mean, mean_err, std, std_err]\n\n")
    f.write("kinematic_variable, {}\n\n".format(kinem))
    f.close()

plt.style.use('grid_multiple_plots')

total_entries = 0

# Determine grid size per page.
n_pdfs = len(eta_ls)-1
n_pages = len(pT_ls) - 1
n_plots_per_page = len(qd0_ls) - 1
rows, cols = ncolsrows_from_nplots(n_plots_per_page)

msg = "Making {} PDF.".format(n_pdfs)
if n_pdfs > 1:
    msg = msg.replace("PDF", "PDFs")
print(msg)
print("  Each PDF will contain {} pages.".format(n_pages))
print("  Each page will contain {} plots ({} x {} grid).\n".format(n_plots_per_page, rows, cols))
print("|eta| regions: {}\n".format(np.round(eta_ls, decimals=2)))
print("pT regions: {}\n".format(np.round(pT_ls, decimals=2)))
print("q*d0 regions: {}".format(np.round(qd0_ls, decimals=4)))
#----------------#
#----- Main -----#
#----------------#
# Loop over eta regions.
for k in range(len(eta_ls)-1):
    eta_min = eta_ls[k]
    eta_max = eta_ls[k+1]
    eta_range = [eta_min, eta_max]

    # For each eta range, make a pdf. 
    extra   = "__{}_eta_{}".format(eta_min, eta_max)
    extra  += "__{:.1f}_pT_{:.1f}_GeV".format(min(pT_ls), max(pT_ls))
    extra = make_str_title_friendly(extra)
    extra += ".pdf"

    fullpath = os.path.join(outdir, filename_base + extra)
    check_overwrite(fullpath, overwrite=overwrite)    

    print("Making PDF {}/{}...".format(k+1, n_pdfs))
    t_start_page = time.perf_counter()
    with PdfPages(fullpath) as pdf:
        # Within 1 eta region, scan the pT regions. 
        for j in range(len(pT_ls)-1):
            pT_min = pT_ls[j]
            pT_max = pT_ls[j+1]
            pT_range = [pT_min, pT_max]

            # Within 1 pT region, scan the q*d0 regions and make 1 page of PDF. 
            t_start = time.perf_counter()
            f = plt.figure()
            for count in range(len(qd0_ls)-1):
                ax = plt.subplot(rows,cols,count+1)
                qd0_min = qd0_ls[count]
                qd0_max = qd0_ls[count+1]
                qd0_range = [qd0_min, qd0_max]

                x_bins, bin_width = make_binning_array(x_bin_info)

                x_label = label_LaTeX_dict[kinem]["independent_label"]
                x_units = label_LaTeX_dict[kinem]["units"]
                y_label = hist_y_label(bin_width, x_units)
                if len(x_units) > 0:
                    x_label += " [{}]".format(x_units)

                cuts  = r"$%.2f < \left| \eta^{\mathrm{REC}} \right| < %.2f$," % (eta_min, eta_max) + "\n"
                cuts += r"$%.1f <$ %s $< %.1f$ GeV," % (pT_min, p_str_latex, pT_max) + "\n"
                cuts += r"$%.4f < q(\mu)*d_{0}^{\mathrm{%s}} < %.4f$, " % (qd0_min, "BS", qd0_max) + "\n"
                cuts += r"$\Delta R < %.3f$, " % (dR_max) + "\n"
                cuts += r"$%.1f < m_{Z} < %.1f$ GeV,  " % (massZ_min_DY, massZ_max_DY) + "\n"
                cuts += r"$%.1f < m_{J/\psi} < %.1f$ GeV" % (massZ_min_Jpsi, massZ_max_Jpsi)

                all_masks_DY = vaex_apply_masks(  vdf_concat_MC_2017_DY,   eta_range, pT_range, qd0_range, massZ_minmax_DY,   dR_max)
                all_masks_Jpsi = vaex_apply_masks(vdf_concat_MC_2017_Jpsi, eta_range, pT_range, qd0_range, massZ_minmax_Jpsi, dR_max)

                total_entries += all_masks_DY.sum() + all_masks_Jpsi.sum()
                ax, bin_vals, bin_edges, stats = make_1D_dist(ax=ax, 
                                                           data=np.append(vdf_concat_MC_2017_DY.evaluate(kinem, selection=all_masks_DY), 
                                                                         vdf_concat_MC_2017_Jpsi.evaluate(kinem, selection=all_masks_Jpsi)
                                                                         ),
                                                              x_limits=x_zoom_range,
                                                              x_bins=x_bins, 
                                                             x_label=x_label, 
                                                             y_label=y_label,
                                                             title="",
                                                             y_max=-1,
                                                            log_scale=False)
                ax.text(0.025, 0.78, cuts, horizontalalignment='left', verticalalignment='center', 
                transform=ax.transAxes, bbox=dict(boxstyle='square', facecolor='white', alpha=0.9))

                if (iter_gaus[0]):
                    # Do iterative fitting procedure.
                    # Produce a dictionary of stats for these fits.
                    # The Gaus coeff is usually around the max height of Gaus curve.
                    best_guess_coeff = np.max(bin_vals)  
                    best_guess_mean = stats[1]
                    best_guess_stdev = stats[3]
                    fit_stats_dict, ax = iterative_fit_gaus(iter_gaus[1], bin_edges, bin_vals, 
                                                        param_guess=[best_guess_coeff, 
                                                                     best_guess_mean, 
                                                                     best_guess_stdev],
                                                        param_bounds=([0,-1000,-100], [999999999,1000,100]),
                                                        ax=ax, draw_on_axes=True, verbose=verbose)
                    # Use plotted kinem as the key for this dict of stats. 

                with open(stats_fullpath, "a") as f:
                    f.write("range_eta, {}\n".format(eta_range))
                    f.write("range_pT, {}\n".format(pT_range))
                    f.write("range_qd0, {}\n".format(qd0_range))
                    f.write("hist_stats, {}\n".format(stats))
                    f.write("fit_stats:\n")
                    for key,val in fit_stats_dict.items():
                        f.write("{}, {}\n".format(key, val))
                    f.write("\n")
                    f.close()

            # End qd0 loop.
            t_end = time.perf_counter()
            plt.tight_layout()
            pdf.savefig()
            plt.close("all")

            msg = "  Page {}/{} made. Time taken: {:.2f} s".format(j+1, n_pages, t_end - t_start)
            print_header_message(msg)
            print()

        # End pT loop. Make next page.
    t_end_page = time.perf_counter()
    print("PDF {}/{} made at:\n  {}".format(k+1, n_pdfs, fullpath))
    print("(in {} s)".format(t_end_page - t_start_page))
    # Save this 1 eta reg, 1 pT reg, and all q*d0 plots on one page. 

# Go to next pT reg and next page.
print("All PDFs created.")

total_muons_original = vdf_concat_MC_2017_DY.count() + vdf_concat_MC_2017_Jpsi.count()
print("Total muons expected: {}".format(total_muons_original))
perdif = perc_diff(total_entries, total_muons_original)
print("Total muons found, inclusive all regions: {} (perc. diff. = {:.2f}%)".format(total_entries, perdif))