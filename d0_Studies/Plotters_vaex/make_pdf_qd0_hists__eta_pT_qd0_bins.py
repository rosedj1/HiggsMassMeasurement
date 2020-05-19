# PURPOSE: Make PDFs of q*d0BS distributions, in given eta bins and pT bins.
# NOTES:   User should check User parameters.
#          One PDF is made per eta bin. 
#          Within each eta bin, all pT bins are analyzed. 
#          Within each pT bin, 8 q*d0BS distributions are made (per page).
#          FIXME: Currently, there must be 9 elements in qd0_ls (8 plots).
# SYNTAX:  python <script>.py
# AUTHOR:  Jake Rosenzweig
# DATE:    2020-05-10

import os
import sys
sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')

import matplotlib.pyplot as plt

from d0_Studies.Plotters.vaex_read_MC_2016_DY_dataframe import vdf_concat
from d0_Utils.d0_dicts import label_LaTeX_dict
from d0_Utils.d0_fns import make_binning_array

from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from PyUtils.Utils_Plotting import hist_y_label, make_1D_dist

from matplotlib.backends.backend_pdf import PdfPages

#---------------------------#
#----- User Parameters -----#
#---------------------------#
outdir = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/qd0_hists/"
filename_base = "Test05_MC_2016_DY_fullstats"
overwrite = False

# Binning.
eta_ls = [0.0, 0.2, 0.4]
pT_ls = [7, 9.8, 13.72, 19.2]
qd0_ls = [-0.003, -0.0015, -0.001, -0.0005, 0.0000, 0.0005, 0.001, 0.0015, 0.003]

qd0_bin_width = 0.00002

dR_cut = 0.008
massZ_ls = [60, 120]

p_str_latex = "$p_{T}$"
#----------------------#
#----- Automatons -----#
#----------------------#
massZ_min = massZ_ls[0]
massZ_max = massZ_ls[1]
x_range = [min(qd0_ls), max(qd0_ls)]

makeDirs(outdir)
plt.style.use('grid_multiple_plots')
#----------------#
#----- Main -----#
#----------------#
for k in range(len(eta_ls)-1):
    eta_min = eta_ls[k]
    eta_max = eta_ls[k+1]

    # For each eta range, make a pdf. 
    extra   = "__{:.1f}_eta_{:.1f}".format(eta_min, eta_max)
    extra  += "__{:.1f}_pT_{:.1f}_GeV".format(min(pT_ls), max(pT_ls))
    extra = make_str_title_friendly(extra)
    extra += ".pdf"

    fullpath = os.path.join(outdir, filename_base + extra)
    check_overwrite(fullpath, overwrite=overwrite)    

    with PdfPages(fullpath) as pdf:
        # Within 1 eta region, scan the pT regions. 
        for j in range(len(pT_ls)-1):
            pT_min = pT_ls[j]
            pT_max = pT_ls[j+1]

            f = plt.figure()
            # Within 1 pT region, scan the q*d0 regions. 
            for count in range(len(qd0_ls)-1):
                ax = plt.subplot(4,2,count+1)
                qd0_min = qd0_ls[count]
                qd0_max = qd0_ls[count+1]

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
            # End qd0 loop.
            plt.tight_layout()
            pdf.savefig()
            plt.close("all")
        # End pT loop.
        print("[INFO] PDF made: {}".format(fullpath))
        # Save this 1 eta reg, 1 pT reg, and all q*d0 plots on one page. 

    # Go to next pT reg and next page.