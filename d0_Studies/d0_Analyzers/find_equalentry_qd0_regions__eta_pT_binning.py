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
import pickle
sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')

import numpy as np

from vaex_Utils.vaex_dataframes import (vdf_MC_2017_DY, vdf_MC_2017_Jpsi, vdf_MC_2017_DY,
                                        prepare_vaex_df, vaex_apply_masks)
from d0_Utils.d0_fns import find_equal_hist_regions_unbinned
from d0_Studies.kinematic_bins import (equal_entry_bin_edges_eta_mod1,
                                        equal_entry_bin_edges_eta_mod2,
                                        equal_entry_bin_edges_eta_mod1_wholenum,
                                        equal_entry_bin_edges_pT_sevenfifths_to1000GeV,
                                        bin_edges_pT_sevenfifths_to1000GeV_wholenum,
                                        )
from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite

#---------------------------#
#----- User Parameters -----#
#---------------------------#
# Samples:
vdf_concat_MC_2017_DY = prepare_vaex_df(vdf_MC_2017_DY)
vdf_concat_MC_2017_Jpsi = prepare_vaex_df(vdf_MC_2017_Jpsi)

outdir = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/find_best_binning__eta_pT_qd0/"
filename_base = "Everythingseemstobeworking_Test06"
overwrite = True

# Binning.
eta_ls = equal_entry_bin_edges_eta_mod1_wholenum[5:7]
pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum[4:7]
qd0_limits = [-0.015, 0.015]

r = 10  # Number of regions to split each q*d0 region into. 
algo = ("at_least", 1000)
round_to_n_decimals = 5
verbose = True

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

total_entries = 0

print("|eta| regions:\n  {}".format(np.round(eta_ls, decimals=2)))
print("pT regions:\n  {}\n".format(np.round(pT_ls, decimals=2)))
#----------------#
#----- Main -----#
#----------------#
# Make only 1 PDF. Each page is a different eta cut.
extra   = "__{:.1f}_eta_{:.1f}".format(min(eta_ls), max(eta_ls))
extra  += "__{:.1f}_pT_{:.1f}_GeV".format(min(pT_ls), max(pT_ls))
extra = make_str_title_friendly(extra)
extra += ".csv"

fullpath = os.path.join(outdir, filename_base + extra)
fullpath_pickle = fullpath.replace(".csv", ".pkl")

check_overwrite(fullpath, overwrite=overwrite)
check_overwrite(fullpath_pickle, overwrite=overwrite)

equal_entry_binedge_dict = {}

with open(fullpath, "w") as myfile:
    # Loop over eta regions.
    for k in range(len(eta_ls)-1):
        eta_min = eta_ls[k]
        eta_max = eta_ls[k+1]
        eta_range = [eta_min, eta_max]

        eta_key_str = "eta_bin_left_edge={}".format(eta_min)
        equal_entry_binedge_dict[eta_key_str] = {}

        # Column names.
        col_str  = "{:<10}, {:<11}, {:<13}, ".format("bins_found","bins_wanted", "muons_per_reg")
        col_str += "\t{}\t\t{}\t\t\t{}\n".format("eta_range,", "pT_range,", "q*d0_bins")
        myfile.write(col_str)

        # Within this eta region, scan the pT regions. 
        for count in range(len(pT_ls)-1):
            pT_min = pT_ls[count]
            pT_max = pT_ls[count+1]
            pT_range = [pT_min, pT_max]

            pT_key_str = "pT_bin_left_edge={}".format(pT_min)

            all_masks_DY = vaex_apply_masks(  vdf_concat_MC_2017_DY,   eta_range, pT_range, qd0_limits, massZ_minmax_DY,   dR_max)
            all_masks_Jpsi = vaex_apply_masks(vdf_concat_MC_2017_Jpsi, eta_range, pT_range, qd0_limits, massZ_minmax_Jpsi, dR_max)

            data = np.append(vdf_concat_MC_2017_DY.evaluate("qd0BS", selection=all_masks_DY), 
                             vdf_concat_MC_2017_Jpsi.evaluate("qd0BS", selection=all_masks_Jpsi))

            n_muons = len(data)
            equalentry_binedge_ls, r_updated = find_equal_hist_regions_unbinned(
                                            data, 
                                            r=r, 
                                            verbose=verbose,
                                            algo=algo,
                                            round_to_n_decimals=round_to_n_decimals)

            # Append to dict.
            equal_entry_binedge_dict[eta_key_str][pT_key_str] = equalentry_binedge_ls

            info  = "{:<10}, {:<11}, {:<13}, ".format(r_updated, r, n_muons // r_updated)
            info += "\t{},\t{},\t\t{}\n".format(eta_range, pT_range, equalentry_binedge_ls)
            myfile.write(info)

            total_entries += n_muons
        # End pT loop. Go to next eta range
        myfile.write('\n')
        print("Finished eta_range = {}".format(eta_range))
    # End eta loop.

with open(fullpath_pickle,'wb') as output:
    pickle.dump(equal_entry_binedge_dict, output, pickle.HIGHEST_PROTOCOL)

total_muons_original = vdf_concat_MC_2017_DY.count() + vdf_concat_MC_2017_Jpsi.count()
print("Total muons expected: {}".format(total_muons_original))
print("Total muons found, inclusive all regions: {}".format(total_entries))