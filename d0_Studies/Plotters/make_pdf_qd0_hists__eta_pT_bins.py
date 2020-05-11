# PURPOSE: Make single PDF of q*d0BS distributions, with no cuts on q*d0BS,
#          in given eta bins and pT bins. Each page of PDF is for a single
#          eta region. Then 4 plots per page are binned in pT. 
# NOTES:   User should check User parameters.
#          One PDF is made in total. 
#          Within each eta bin, all pT bins are analyzed. 
#          A single q*d0 hist is made for each pT bin.
#          Each page looks best when there are 8 plots per page, 
#          i.e. 9 elements in pT_ls. 
# SYNTAX:  python <script>.py
# AUTHOR:  Jake Rosenzweig
# DATE:    2020-05-11

import os
import sys
import math
sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')

import matplotlib.pyplot as plt

from d0_Studies.Plotters.vaex_read_MC_2016_DY_dataframe import vdf_concat
from d0_Studies.kinematic_bins import (equal_entry_bin_edges_pT,
                                       equal_entry_bin_edges_pT_mod2,
                                       equal_entry_bin_edges_pT_sevenfifths,
                                       equal_entry_bin_edges_pT_sevenfifths_mod,
                                       equal_entry_bin_edges_eta, 
                                       equal_entry_bin_edges_eta_mod)
from d0_Utils.d0_dicts import label_LaTeX_dict
from d0_Utils.d0_fns import make_binning_array

from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from PyUtils.Utils_Plotting import hist_y_label, make_1D_dist, ncolsrows_from_nplots

from matplotlib.backends.backend_pdf import PdfPages

#---------------------------#
#----- User Parameters -----#
#---------------------------#
outdir = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/qd0_hists/"
filename_base = "MC_2016_DY_fullstats_no_qd0cuts_modpTbins_TEST09"
overwrite = False

# Binning.
eta_ls = equal_entry_bin_edges_eta_mod
pT_ls = equal_entry_bin_edges_pT_mod2
qd0_ls = [-0.020, 0.020]
qd0_bin_width = 0.0002

x_range = [-0.022, 0.022]#[min(qd0_ls), max(qd0_ls)]

dR_cut = 0.008
massZ_ls = [60, 120]

p_str_latex = "$p_{T}$"
#----------------------#
#----- Automatons -----#
#----------------------#
massZ_min = massZ_ls[0]
massZ_max = massZ_ls[1]

makeDirs(outdir)
plt.style.use('grid_multiple_plots')

qd0_min = qd0_ls[0]
qd0_max = qd0_ls[1]

n_plots = len(pT_ls) - 1
rows, cols = ncolsrows_from_nplots(n_plots)
print("[INFO] Making a {} page PDF.".format(len(eta_ls)-1))
print("[INFO] Making {} plots ({} x {} grid per page)".format(n_plots, rows, cols))
#----------------#
#----- Main -----#
#----------------#
# Make only 1 PDF. Each page is a different eta cut.
extra   = "__{:.1f}_eta_{:.1f}".format(min(eta_ls), min(eta_ls))
extra  += "__{:.1f}_pT_{:.1f}_GeV".format(min(pT_ls), max(pT_ls))
extra = make_str_title_friendly(extra)
extra += ".pdf"

fullpath = os.path.join(outdir, filename_base + extra)
check_overwrite(fullpath, overwrite=overwrite)    
with PdfPages(fullpath) as pdf:
    for k in range(len(eta_ls)-1):
        eta_min = eta_ls[k]
        eta_max = eta_ls[k+1]

        f = plt.figure()
        # Within this eta region, scan the pT regions. 
        for count in range(len(pT_ls)-1):
            ax = plt.subplot(rows,cols,count+1)
            pT_min = pT_ls[count]
            pT_max = pT_ls[count+1]

            # Selections.
            mask_eta = (eta_min < vdf_concat["eta"]) & (vdf_concat["eta"] < eta_max)
            mask_pT = (pT_min < vdf_concat["pT"]) & (vdf_concat["pT"] < pT_max)
            mask_qd0 = (qd0_min < vdf_concat["qd0BS"]) & (vdf_concat["qd0BS"] < qd0_max)
            mask_massZ = (60 < vdf_concat["massZ"]) & (vdf_concat["massZ"] < 120)
            mask_dR = (vdf_concat["delta_R"] < 0.008)

            all_masks = mask_pT & mask_eta & mask_qd0 & mask_massZ & mask_dR

            x_bin_limits = [qd0_min, qd0_max, qd0_bin_width]

            x_bins, binwidth = make_binning_array(x_bin_limits)

            x_label = label_LaTeX_dict['qd0BS1']["independent_label"]
            x_units = label_LaTeX_dict['qd0BS1']["units"]
            y_label = hist_y_label(binwidth, x_units)
            if len(x_units) > 0:
                x_label += " [{}]".format(x_units)

            cuts  = r"$%.2f < \left| \eta^{\mathrm{REC}} \right| < %.2f$," % (eta_min, eta_max) + "\n"
            cuts += r"$%.1f <$ %s $< %.1f$ GeV," % (pT_min, p_str_latex, pT_max) + "\n"
            cuts += r"$%.4f < q(\mu)*d_{0}^{\mathrm{%s}} < %.4f$, " % (qd0_min, "BS", qd0_max) + "\n"
            cuts += r"$\Delta R < %.3f$, " % (dR_cut) + "\n"
            cuts += r"$%.1f < m_{\mu\mu} < %.1f$ GeV" % (massZ_min, massZ_max)

            ax, bin_vals, bin_edges, stats = make_1D_dist(ax=ax, 
                                                          data=vdf_concat.evaluate("qd0BS", selection=all_masks), 
                                                          x_limits=x_range,
                                                          x_bins=x_bins, 
                                                         x_label=x_label, 
                                                         y_label=y_label,
                                                         title="",
                                                         y_max=-1,
                                                        log_scale=False)
            ax.text(0.025, 0.78, cuts, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes,
                   bbox=dict(boxstyle='square', facecolor='white', alpha=0.9))
        # End pT loop.
        plt.tight_layout()
        pdf.savefig()
        plt.close("all")
        print("[INFO] Page {} made.".format(k+1))
    # End eta loop.
    print("[INFO] PDF made at:\n  {}".format(fullpath))