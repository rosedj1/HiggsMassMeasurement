"""
# Purpose: 
#   Find the bin edges of q*d0 distributions such that all regions 
#   of the dist. have equal entries.
#   Each q*d0 distribution will have eta and pT cuts applied. 
#   User specifies the number of equal-entry regions to split each dist up into.
#   Bin edges of equal-entry regions are recorded in a .csv file, 
#   and a .pkl file. A dictionary is saved to the .pkl. 
# Syntax:  
#   python script.py 
#   python script.py > output.txt
# Notes:   
#   This code runs on DY, J/psi, and DY+J/psi, treating each individually.
#   User puts the number of regions (r) to split each q*d0 distribution, 
#   such that each region will have equal entries in it.
#   User can also request that "at_least" N entries are found
#   per region. If r did not give N, then r is decremented 
#   until either N is found or r == 2. 
#   Make sure to check all the parameters in "User Parameters".
#   Should be used with Python 3.X.
# Author:  Jake Rosenzweig
# Updated: 2020-05-24
"""
import os
import sys
import pickle

import numpy as np

from vaex_Utils.vaex_dataframes import (vdf_MC_2017_DY, vdf_MC_2017_Jpsi, vdf_MC_2017_DY,
                                        prepare_vaex_df, vaex_apply_masks)
from d0_Utils.d0_fns import find_equal_hist_regions_unbinned
from d0_Studies.kinematic_bins import (equal_entry_bin_edges_eta_mod1_wholenum,
                                        bin_edges_pT_sevenfifths_to1000GeV_wholenum)
from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite

#---------------------------#
#----- User Parameters -----#
#---------------------------#
# Samples:
vdf_concat_MC_2017_DY = prepare_vaex_df(vdf_MC_2017_DY)
vdf_concat_MC_2017_Jpsi = prepare_vaex_df(vdf_MC_2017_Jpsi)

# outdir = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/find_best_binning__eta_pT_qd0/"
# Where to save pkl and csv files.
outdir = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/"

filename_base = "fullscan_individ_and_comb_samples"
write_to_pickle = True
overwrite = False
verbose = True

# Binning.
eta_ls = equal_entry_bin_edges_eta_mod1_wholenum
pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum
# eta_ls = [1.4, 1.6]
# pT_ls = [30.0, 40.0, 50.0]
qd0_limits = [-0.015, 0.015]

r = 12  # Number of regions to split each q*d0 region into. 
algo = ("at_least", 2000)

dR_max = 0.008
massZ_minmax_DY = [60, 120]
massZ_minmax_Jpsi = [2.9, 3.3]

#----------------------#
#----- Automatons -----#
#----------------------#
massZ_min_DY = massZ_minmax_DY[0]
massZ_max_DY = massZ_minmax_DY[1]
massZ_min_Jpsi = massZ_minmax_Jpsi[0]
massZ_max_Jpsi = massZ_minmax_Jpsi[1]

makeDirs(outdir)

total_entries = 0

if (verbose):
    print(f"|eta| regions:\n  {np.round(eta_ls, decimals=2)}")
    print(f"pT regions:\n  {np.round(pT_ls, decimals=2)}\n")
#----------------#
#----- Main -----#
#----------------#
extra = (
    f"_{r}reg"
    f"with{algo[1]}perreg"
    f"__{min(eta_ls):.1f}_eta_{max(eta_ls):.1f}"
    f"__{min(pT_ls):.1f}_pT_{max(pT_ls):.1f}_GeV"
)
extra = make_str_title_friendly(extra)
extra += ".csv"

fullpath_csv = os.path.join(outdir, filename_base + extra)
fullpath_pickle = fullpath_csv.replace(".csv", ".pkl")

check_overwrite(fullpath_csv, overwrite=overwrite)
check_overwrite(fullpath_pickle, overwrite=overwrite)

equal_entry_binedge_dict = {}

with open(fullpath_csv, "w") as myfile:
    myfile.write(f"all_eta_bins : {eta_ls}\n")
    myfile.write(f"all_pT_bins  : {pT_ls}\n\n")

    # Loop over eta regions.
    for k in range(len(eta_ls)-1):
        eta_min = eta_ls[k]
        eta_max = eta_ls[k+1]
        eta_range = [eta_min, eta_max]

        eta_key_str = f"eta_bin_left_edge={eta_min}"
        equal_entry_binedge_dict[eta_key_str] = {}

        # Column names.
        col_str  = (
            f"{'sample':<7}, {'bins_found':<2}, {'bins_wanted':<2}, {'muons_per_reg':<10}, "
            f"\teta_range,\t\tpT_range,\t\tq*d0_bins\n"
        )
        myfile.write(col_str)

        # Within this eta region, scan the pT regions. 
        max_count = len(pT_ls)-1
        for count in range(max_count):
            pT_min = pT_ls[count]
            pT_max = pT_ls[count+1]
            pT_range = [pT_min, pT_max]

            pT_key_str = "pT_bin_left_edge={}".format(pT_min)
            equal_entry_binedge_dict[eta_key_str][pT_key_str] = {}
            # No restriction on q*d0.
            all_masks_DY = vaex_apply_masks(  vdf_concat_MC_2017_DY,   eta_range, pT_range, qd0_limits, massZ_minmax_DY,   dR_max)
            all_masks_Jpsi = vaex_apply_masks(vdf_concat_MC_2017_Jpsi, eta_range, pT_range, qd0_limits, massZ_minmax_Jpsi, dR_max)

            n_tot_DY = all_masks_DY.count()
            n_tot_Jpsi = all_masks_Jpsi.count()
            n_sel_DY = all_masks_DY.sum()
            n_sel_Jpsi = all_masks_Jpsi.sum()
            # A couple of checks.
            assert vdf_concat_MC_2017_DY.count() == n_tot_DY
            assert vdf_concat_MC_2017_Jpsi.count() == n_tot_Jpsi

            if (verbose):
                print(f"pT_loop_count: {count}/{max_count}")
                print(f"  eta_range={eta_range}")
                print(f"  pT_range={pT_range}")

                print(f"Total DY events: {n_tot_DY}")
                print(f"Total Jpsi events: {n_tot_Jpsi}")
                print(f"DY events after selection: {n_sel_DY} (eff. = {n_sel_DY/float(n_tot_DY)*100.:.4f})%")
                print(f"Jpsi events after selection: {n_sel_Jpsi} (eff. = {n_sel_Jpsi/float(n_tot_Jpsi)*100.:.4f})%")

            # Get ALL q*d0 values.
            qd0_arr_DY = vdf_concat_MC_2017_DY.evaluate("qd0BS")  # These are numpy arrays.
            qd0_arr_Jpsi = vdf_concat_MC_2017_Jpsi.evaluate("qd0BS")

            # Remember this is within just one (eta, pT) region.
            qd0_arr_sel_DY = qd0_arr_DY[all_masks_DY.values]
            qd0_arr_sel_Jpsi = qd0_arr_Jpsi[all_masks_Jpsi.values]
            qd0_arr_sel_comb = np.concatenate((qd0_arr_sel_DY, qd0_arr_sel_Jpsi))

            # Find equal-entry q*d0 regions for individual samples and for combined samples.
            equalentry_binedge_ls_DY, r_updated_DY = find_equal_hist_regions_unbinned(
                                            qd0_arr_sel_DY, 
                                            r=r, 
                                            verbose=verbose,
                                            algo=algo)

            equalentry_binedge_ls_Jpsi, r_updated_Jpsi = find_equal_hist_regions_unbinned(
                                            qd0_arr_sel_Jpsi, 
                                            r=r, 
                                            verbose=verbose,
                                            algo=algo)

            equalentry_binedge_ls_comb, r_updated_comb = find_equal_hist_regions_unbinned(
                                            qd0_arr_sel_comb, 
                                            r=r, 
                                            verbose=verbose,
                                            algo=algo)

            # Append to dict.
            equal_entry_binedge_dict[eta_key_str][pT_key_str]["equalentry_qd0ls_DY"] = equalentry_binedge_ls_DY
            equal_entry_binedge_dict[eta_key_str][pT_key_str]["equalentry_qd0ls_Jpsi"] = equalentry_binedge_ls_Jpsi
            equal_entry_binedge_dict[eta_key_str][pT_key_str]["equalentry_qd0ls_DY+Jpsi"] = equalentry_binedge_ls_comb

            n_muons_DY = len(qd0_arr_sel_DY)
            n_muons_Jpsi = len(qd0_arr_sel_Jpsi)
            n_muons_comb = len(qd0_arr_sel_comb)
            info  = (
                f"{'DY':<7}, {r_updated_DY:<2}, {r:<2}, {n_muons_DY // r_updated_DY:<10}, "
                f"\t{eta_range},\t\t{pT_range},\t\t{equalentry_binedge_ls_DY}\n"
                f"{'Jpsi':<7}, {r_updated_Jpsi:<2}, {r:<2}, {n_muons_Jpsi // r_updated_Jpsi:<10}, "
                f"\t{eta_range},\t\t{pT_range},\t\t{equalentry_binedge_ls_Jpsi}\n"
                f"{'DY+Jpsi':<7}, {r_updated_comb:<2}, {r:<2}, {n_muons_comb // r_updated_comb:<10}, "
                f"\t{eta_range},\t\t{pT_range},\t\t{equalentry_binedge_ls_comb}\n"
            )
            myfile.write(info)

            total_entries += n_muons_DY + n_muons_Jpsi
            print(f"total_entries up to this point: {total_entries}")
        # End pT loop. Go to next eta range
        myfile.write('\n')
        print(f"Finished eta_range = {eta_range}\n\n")
    # End eta loop.
    
print(f"[INFO] q*d0 bin edge info written to csv file:\n{fullpath_csv}")

if (write_to_pickle):
    with open(fullpath_pickle,'wb') as output:
        pickle.dump(equal_entry_binedge_dict, output, pickle.HIGHEST_PROTOCOL)
    print(f"[INFO] eta, pT, q*d0 bin dict written to pickle file:\n{fullpath_pickle}\n")

total_muons_original = vdf_concat_MC_2017_DY.count() + vdf_concat_MC_2017_Jpsi.count()
perc = total_entries / float(total_muons_original) * 100.
print(
    f"Total muons found before selections: {total_muons_original}\n"
    f"Total muons found after selections (inclusive all regions): {total_entries} ({perc:.2f}%)"
)