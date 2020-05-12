import time
import os
import sys
import math
sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')

import matplotlib.pyplot as plt

from vaex_Utils.vaex_dataframes import (vdf_MC_2017_DY, vdf_MC_2017_Jpsi, vdf_MC_2017_DY,
                                        prepare_vaex_df, vaex_apply_masks)
from d0_Studies.kinematic_bins import (equal_entry_bin_edges_pT,
                                       equal_entry_bin_edges_pT_mod2,
                                       equal_entry_bin_edges_pT_sevenfifths,
                                       equal_entry_bin_edges_pT_sevenfifths_mod,
                                       equal_entry_bin_edges_eta, 
                                       equal_entry_bin_edges_eta_mod1, 
                                       equal_entry_bin_edges_eta_mod2,)
from d0_Utils.d0_dicts import label_LaTeX_dict
from d0_Utils.d0_fns import make_binning_array

from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from PyUtils.Utils_Plotting import hist_y_label, make_1D_dist, ncolsrows_from_nplots

from matplotlib.backends.backend_pdf import PdfPages

#---------------------------#
#----- User Parameters -----#
#---------------------------#
# Samples:
vdf_concat_MC_2017_DY = prepare_vaex_df(vdf_MC_2017_DY)

outdir = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/qd0_hists/"
filename_base = "MC2017DY_fullstats_sevenfifths_pT__no_qd0cuts_TESTdeploy02"
overwrite = False

# Binning.
eta_ls = equal_entry_bin_edges_eta_mod2[:4]#[0.0, 0.2]#
pT_ls = equal_entry_bin_edges_pT_sevenfifths#equal_entry_bin_edges_pT_mod2
qd0_bin_info = [-0.02, 0.02, 0.0002]

x_range = [-0.022, 0.022]

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

makeDirs(outdir)
plt.style.use('grid_multiple_plots')

qd0_min = qd0_bin_info[0]
qd0_max = qd0_bin_info[1]
qd0_range = qd0_bin_info[0:2]
x_bins, bin_width = make_binning_array(qd0_bin_info)
    
n_plots = len(pT_ls) - 1
rows, cols = ncolsrows_from_nplots(n_plots)
print("[INFO] Making a {}-page PDF.".format(len(eta_ls)-1))
print("[INFO] Making {} plots ({} x {} grid per page)".format(n_plots, rows, cols))
#----------------#
#----- Main -----#
#----------------#
# Make only 1 PDF. Each page is a different eta cut.
extra   = "__{:.1f}_eta_{:.1f}".format(min(eta_ls), max(eta_ls))
extra  += "__{:.1f}_pT_{:.1f}_GeV".format(min(pT_ls), max(pT_ls))
extra = make_str_title_friendly(extra)
extra += ".pdf"

fullpath = os.path.join(outdir, filename_base + extra)
check_overwrite(fullpath, overwrite=overwrite)    

with PdfPages(fullpath) as pdf:
    # Loop over eta regions. 
    for k in range(len(eta_ls)-1):
        eta_min = eta_ls[k]
        eta_max = eta_ls[k+1]

        eta_range = [eta_min, eta_max]
        
        # Within this eta region, scan the pT regions. 
        t_start = time.perf_counter()
        f = plt.figure()
        for count in range(len(pT_ls)-1):
            ax = plt.subplot(rows,cols,count+1)
            pT_min = pT_ls[count]
            pT_max = pT_ls[count+1]
            pT_range = [pT_min, pT_max]

            x_label = label_LaTeX_dict['qd0BS1']["independent_label"]
            x_units = label_LaTeX_dict['qd0BS1']["units"]
            y_label = hist_y_label(bin_width, x_units)
            if len(x_units) > 0:
                x_label += " [{}]".format(x_units)

            cuts  = r"$%.2f < \left| \eta^{\mathrm{REC}} \right| < %.2f$," % (eta_min, eta_max) + "\n"
            cuts += r"$%.1f <$ %s $< %.1f$ GeV," % (pT_min, p_str_latex, pT_max) + "\n"
            cuts += r"$%.4f < q(\mu)*d_{0}^{\mathrm{%s}} < %.4f$, " % (qd0_min, "BS", qd0_max) + "\n"
            cuts += r"$\Delta R < %.3f$, " % (dR_max) + "\n"
            cuts += r"$%.1f < m_{Z} < %.1f$ GeV,  " % (massZ_min_DY, massZ_max_DY) + "\n"
            cuts += r"$%.1f < m_{J/\psi} < %.1f$ GeV" % (massZ_min_Jpsi, massZ_max_Jpsi)

            all_masks_DY = vaex_apply_masks(vdf_concat_MC_2017_DY, 
                                            eta_range, pT_range, qd0_range, massZ_minmax_DY, 
                                            dR_max)
            ax, bin_vals, bin_edges, stats = make_1D_dist(ax=ax, 
                                                          data=vdf_concat_MC_2017_DY.evaluate("qd0BS", selection=all_masks_DY),
                                                          x_limits=x_range,
                                                          x_bins=x_bins, 
                                                          x_label=x_label, 
                                                          y_label=y_label,
                                                          title="",
                                                          y_max=-1,
                                                          log_scale=False)
            ax.text(0.025, 0.78, cuts, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes,
                          bbox=dict(boxstyle='square', facecolor='white', alpha=0.9))
            t_end = time.perf_counter()
        # End pT loop.
        
        plt.tight_layout()
        pdf.savefig()
        plt.close("all")
        print("[INFO] Page {} made. Time taken: {:.2f} s".format(k+1, t_end - t_start))

    # End eta loop.
    print("[INFO] PDF made at:\n  {}".format(fullpath))