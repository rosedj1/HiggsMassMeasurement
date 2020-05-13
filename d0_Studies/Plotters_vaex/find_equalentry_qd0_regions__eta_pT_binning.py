"""
# Purpose: Write information to a file that describes equal-entry regions in q*d0 distributions.
# Syntax:  python <script.py>
# Notes:   User puts the number of regions (r) to split each q*d0 distribution, 
#          such that each region will have equal entries in it.
#          Make sure to check all the parameters in "User Parameters".
# Author:  Jake Rosenzweig
# Date:    2020-05-13
"""
import os
import sys
sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')

import numpy as np

from vaex_Utils.vaex_dataframes import (vdf_MC_2017_DY, vdf_MC_2017_Jpsi, vdf_MC_2017_DY,
                                        prepare_vaex_df, vaex_apply_masks)
from d0_Studies.kinematic_bins import (equal_entry_bin_edges_eta_mod1,
                                        equal_entry_bin_edges_eta_mod2,
                                        equal_entry_bin_edges_eta_mod1_wholenum,
                                        equal_entry_bin_edges_pT_sevenfifths_to1000GeV,
                                        bin_edges_pT_sevenfifths_to1000GeV_wholenum,
                                        )
from d0_Utils.d0_fns import find_equal_hist_regions_unbinned

from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite

#---------------------------#
#----- User Parameters -----#
#---------------------------#
# Samples:
vdf_concat_MC_2017_DY = prepare_vaex_df(vdf_MC_2017_DY)
vdf_concat_MC_2017_Jpsi = prepare_vaex_df(vdf_MC_2017_Jpsi)

# outdir = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/find_best_binning__eta_pT_qd0/"
outdir = "/Users/Jake/HiggsMassMeasurement/d0_Studies/Plotters_vaex/"
filename_base = "Test04"
overwrite = True

# Binning.
eta_ls = equal_entry_bin_edges_eta_mod1_wholenum[:2]
pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum[:6]
qd0_bin_info = [-0.015, 0.015, 0.0002]

dR_max = 0.008
massZ_minmax_DY = [60, 120]
massZ_minmax_Jpsi = [2.9, 3.3]

r = 10  # Number of regions to split each q*d0 region into. 
round_to_n_decimals = 5
verbose = False

p_str_latex = "$p_{T}$"

#----------------------#
#----- Automatons -----#
#----------------------#
massZ_min_DY = massZ_minmax_DY[0]
massZ_max_DY = massZ_minmax_DY[1]
massZ_min_Jpsi = massZ_minmax_Jpsi[0]
massZ_max_Jpsi = massZ_minmax_Jpsi[1]

makeDirs(outdir)

qd0_min = qd0_bin_info[0]
qd0_max = qd0_bin_info[1]
qd0_range = qd0_bin_info[0:2]

total_entries = 0

print("|eta| regions:\n  {}".format(np.round(eta_ls, decimals=2)))
print("pT regions:\n  {}".format(np.round(pT_ls, decimals=2)))
#----------------#
#----- Main -----#
#----------------#
# Make only 1 PDF. Each page is a different eta cut.
extra   = "__{:.1f}_eta_{:.1f}".format(min(eta_ls), max(eta_ls))
extra  += "__{:.1f}_pT_{:.1f}_GeV".format(min(pT_ls), max(pT_ls))
extra = make_str_title_friendly(extra)
extra += ".csv"

fullpath = os.path.join(outdir, filename_base + extra)
check_overwrite(fullpath, overwrite=overwrite)

with open(fullpath, "w") as myfile:

    # Loop over eta regions.
    for k in range(len(eta_ls)-1):
        eta_min = eta_ls[k]
        eta_max = eta_ls[k+1]
        eta_range = [eta_min, eta_max]

        # Within this eta region, scan the pT regions. 
        for count in range(len(pT_ls)-1):
            pT_min = pT_ls[count]
            pT_max = pT_ls[count+1]
            pT_range = [pT_min, pT_max]

            all_masks_DY = vaex_apply_masks(  vdf_concat_MC_2017_DY,   eta_range, pT_range, qd0_range, massZ_minmax_DY,   dR_max)
            all_masks_Jpsi = vaex_apply_masks(vdf_concat_MC_2017_Jpsi, eta_range, pT_range, qd0_range, massZ_minmax_Jpsi, dR_max)

            data = np.append(vdf_concat_MC_2017_DY.evaluate("qd0BS", selection=all_masks_DY), 
                             vdf_concat_MC_2017_Jpsi.evaluate("qd0BS", selection=all_masks_Jpsi))

            n_muons = len(data)
            equalentry_binedge_ls = find_equal_hist_regions_unbinned(data, r=r, verbose=verbose, round_to_n_decimals=round_to_n_decimals)

            info  = "n_per_reg, {},\t".format(n_muons // r)
            info += "eta, {},\t".format(eta_range)
            info += "pT, {},\t\t".format(pT_range)
            info += "q*d0_bins, {}\n".format(equalentry_binedge_ls)
            myfile.write(info)

            total_entries += n_muons
        # End pT loop. Go to next eta range
        myfile.write('\n')
        print("Finished eta_range = {}".format(eta_range))
    # End eta loop.

total_muons_original = vdf_concat_MC_2017_DY.count() + vdf_concat_MC_2017_Jpsi.count()
print("Total muons expected: {}".format(total_muons_original))
print("Total muons found, inclusive all regions: {}".format(total_entries))