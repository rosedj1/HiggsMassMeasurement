"""
# PURPOSE: 
#   Save the delta_pT/pT data for each (eta, pT, q*d0) cube. 
#   Also saves the statistics of pT and q*d0 (like the avg, std, and errors)
#   into a dict, which is stored in a '.pkl'.
#   All these data are to be eventually fed into a RooFit unbinned fit. 
# NOTES:   
#   User must specify the input pkl file which contains the equal-entry bin edges:
#     pickled_dict = {eta_bin_left_edge : {
                        pT_bin_left_edge : {
                            "equalentry_qd0ls_DY" : [],
                            "equalentry_qd0ls_Jpsi" : [],
                            "equalentry_qd0ls_DY+Jpsi" : [],
                        }
                    }
#   User should check User Parameters.
# SYNTAX:  python script.py
# AUTHOR:  Jake Rosenzweig
# UPDATED: 2020-05-25
"""
import os
import sys
import pickle
import argparse 

import numpy as np

# Local imports.
from vaex_Utils.vaex_dataframes import (vdf_MC_2017_DY, vdf_MC_2017_Jpsi, vdf_MC_2016_DY,
                                        prepare_vaex_df, vaex_apply_masks)
from d0_Studies.kinematic_bins import (equal_entry_bin_edges_eta_mod1)
from PyUtils.Utils_Physics import perc_diff
from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from PyUtils.Utils_Plotting import get_stats_1Dhist

def ParseOption():
    parser = argparse.ArgumentParser(description='submit all')
    parser.add_argument('--pklfilename', dest='filename_base', type=str, help='') 
    parser.add_argument('--verbose', dest='verbose', type=int, default=1, help='')
    parser.add_argument('--overwrite', dest='overwrite', type=int, default=0, help='')  
    
    args = parser.parse_args()                                                                                         
    return args          
                                                                                                         
# args = ParseOption()          

#---------------------------#
#----- User Parameters -----#
#---------------------------#
# Samples with muons treated "individually":
vdf_concat_MC_2017_DY = prepare_vaex_df(vdf_MC_2017_DY)
vdf_concat_MC_2017_Jpsi = prepare_vaex_df(vdf_MC_2017_Jpsi)

# Dictionary which contains equal-entry q*d0 bin edges.
# inpath_3Dbins_pickle_dict = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/fullscan_individ_and_comb_samples_12regwith2000perreg__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"
inpath_3Dbins_pickle_dict = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/20200525_fullstats_12regwith2000perreg__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"

#-----Below is a big boy. -----#
# inpath_3Dbins_pickle_dict = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/find_best_binning__eta_pT_qd0/equalentry_qd0_binedges__20_regions_max__atleast1000entriesperregion__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"
# outdir_kinem_arrays_dir = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/kinbin3D_pkls/"
outdir_kinem_arrays_dir = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/"

# Filename of outdict is automatic.
overwrite = False
verbose = True

# Cuts to make.
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

# Prepare directories and file names. 
makeDirs(outdir_kinem_arrays_dir)

with open(inpath_3Dbins_pickle_dict, "rb") as f:
    equalentry_binedge_dict = pickle.load(f)

# Binning automatically detected from input pickle file.
eta_ls = equalentry_binedge_dict["all_eta_bins"]
pT_ls = equalentry_binedge_dict["all_pT_bins"]

data_dict = {
    "all_eta_bins" : eta_ls,
    "all_pT_bins" : pT_ls,
    "sample_mod_ls" : [],
}

print(
    f"eta_ls: {eta_ls}\n",
    f"pT_ls: {pT_ls}",
)

filename = os.path.split(inpath_3Dbins_pickle_dict)[1]
filename = filename.replace(".pkl", "_datadict.pkl")
fullpath_kinem_arrays_pkl = os.path.join(outdir_kinem_arrays_dir, filename) 
check_overwrite(fullpath_kinem_arrays_pkl, overwrite=overwrite) 
# If file is there and you wish to overwrite, 
# make sure not to append to existing file; just get rid of it. 
if (overwrite) and os.path.exists(fullpath_kinem_arrays_pkl):
    os.remove(fullpath_kinem_arrays_pkl)

total_num_muons_DY = total_num_muons_Jpsi = total_entries = 0

# Unpack kinematic data.
eta_arr_DY     = vdf_concat_MC_2017_DY.evaluate("eta")
eta_arr_Jpsi   = vdf_concat_MC_2017_Jpsi.evaluate("eta")
pT_arr_DY      = vdf_concat_MC_2017_DY.evaluate("pT")
pT_arr_Jpsi    = vdf_concat_MC_2017_Jpsi.evaluate("pT")
qd0_arr_DY     = vdf_concat_MC_2017_DY.evaluate("qd0BS")
qd0_arr_Jpsi   = vdf_concat_MC_2017_Jpsi.evaluate("qd0BS")
massZ_arr_DY   = vdf_concat_MC_2017_DY.evaluate("massZ")
massZ_arr_Jpsi = vdf_concat_MC_2017_Jpsi.evaluate("massZ")
dR_arr_DY      = vdf_concat_MC_2017_DY.evaluate("delta_R")
dR_arr_Jpsi    = vdf_concat_MC_2017_Jpsi.evaluate("delta_R")

# Prep the masks that don't change. 
mask_massZ_DY   = (massZ_min_DY < massZ_arr_DY) & (massZ_arr_DY < massZ_max_DY)
mask_massZ_Jpsi = (massZ_min_Jpsi < massZ_arr_Jpsi) & (massZ_arr_Jpsi < massZ_max_Jpsi)
mask_dR_DY   = (dR_arr_DY < dR_max)
mask_dR_Jpsi = (dR_arr_Jpsi < dR_max)

# Kinematic to be plotted in each KinBin3D.
# Kinematic data to be stored.
# Should not contain 1 or 2. 
# Acceptable values found in prepare_vaex_df().
kinem_arr_DY_dpToverGenpT = vdf_concat_MC_2017_DY.evaluate("delta_pToverGenpT")
kinem_arr_DY_pT = vdf_concat_MC_2017_DY.evaluate("pT")
kinem_arr_DY_qd0 =  vdf_concat_MC_2017_DY.evaluate("qd0BS")
kinem_arr_Jpsi_dpToverGenpT = vdf_concat_MC_2017_Jpsi.evaluate("delta_pToverGenpT")
kinem_arr_Jpsi_pT = vdf_concat_MC_2017_Jpsi.evaluate("pT")
kinem_arr_Jpsi_qd0 =  vdf_concat_MC_2017_Jpsi.evaluate("qd0BS")

#----------------#
#----- Main -----#
#----------------#
# Loop over eta regions.
for k in range(len(eta_ls)-1):
    eta_min = eta_ls[k]
    eta_max = eta_ls[k+1]
    eta_range = [eta_min, eta_max]
    print(f"eta loop k={k}")

    eta_key = "eta_bin_left_edge={}".format(eta_min)
    data_dict[eta_key] = {}

    # Within 1 eta region, scan the pT regions. 
    for j in range(len(pT_ls)-1):
        pT_min = pT_ls[j]
        pT_max = pT_ls[j+1]
        pT_range = [pT_min, pT_max]
        print(f"pT loop j={j}")

        pT_key = "pT_bin_left_edge={}".format(pT_min)
        data_dict[eta_key][pT_key] = {}

        # Get all equal-entry qd0 bin edges and acquire kinem data. 
        # Remember that there are three different samples, each with their
        # own qd0 ls: DY, J/psi, DY+J/psi
        # Also remember that the lengths of the different qd0 ls may differ!
        qd0_dict = equalentry_binedge_dict[eta_key][pT_key]
        
        for sample, qd0_ls in qd0_dict.items():
            sample_mod = sample.strip("equalentry_qd0ls_")
            data_dict[eta_key][pT_key][sample_mod] = {}
            if sample_mod not in data_dict["sample_mod_ls"]:
                data_dict["sample_mod_ls"].append(sample_mod) 
            print("#" + "-"*10 + "#\n")
            print(f"Running over {sample_mod} with qd0_ls:\n{qd0_ls}")

            for count in range(len(qd0_ls)-1):
                qd0_min = qd0_ls[count]
                qd0_max = qd0_ls[count+1]
                qd0_range = [qd0_min, qd0_max]

                qd0_key = "qd0_bin_left_edge={}".format(qd0_min)
                data_dict[eta_key][pT_key][sample_mod][qd0_key] = {}

                if (verbose):
                    print(f"\nqd0_loop={count}")
                    print(f"sample_mod={sample_mod}")
                    print(f"eta_range={eta_range}")
                    print(f"pT_range={pT_range}")
                    print(f"qd0_range={qd0_range}")

                # Make masks, apply cuts, get data.
                try:
                    # Prepare masks.
                    # These must be labelled with "DY" since they are of DY's size.
                    if "DY" in sample_mod:
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

                    if "Jpsi" in sample_mod:
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
                    print(f"[WARNING] sample_mod={sample_mod} threw a TypeError.")
                    print("There were probably 0 muons found in this region. Continuing anyway.")
                    data_dict[eta_key][pT_key][sample_mod][qd0_key]["dpToverpT_vals"] = [None]
                    data_dict[eta_key][pT_key][sample_mod][qd0_key]["stats_ls_pT"] = [None]
                    data_dict[eta_key][pT_key][sample_mod][qd0_key]["stats_ls_qd0"] = [None]
                    continue

                # Collect the data with selections to be stored in dict.
                if sample_mod == "DY":   
                    mask = all_masks_DY
                    data = data_DY_dpToverGenpT
                    stats_ls_pT = get_stats_1Dhist(data_DY_pT)
                    stats_ls_qd0 = get_stats_1Dhist(data_DY_qd0)

                    n_DY = len(data)
                    total_num_muons_DY += n_DY
                    if (verbose): 
                        print(f"  Number of DY muons in this cube: {n_DY}")
                        print(f"  Cumulative number of DY muons:   {total_num_muons_DY}")
                elif sample_mod == "Jpsi": 
                    mask = all_masks_Jpsi
                    data = data_Jpsi_dpToverGenpT
                    stats_ls_pT = get_stats_1Dhist(data_Jpsi_pT)
                    stats_ls_qd0 = get_stats_1Dhist(data_Jpsi_qd0)

                    n_Jpsi = len(data)
                    total_num_muons_Jpsi += n_Jpsi
                    if (verbose): 
                        print(f"  Number of Jpsi muons in this cube: {n_Jpsi}")
                        print(f"  Cumulative number of Jpsi muons:   {total_num_muons_Jpsi}")
                elif sample_mod == "DY+Jpsi":
                    print("Concatenating DY+Jpsi data together")
                    mask = np.concatenate((all_masks_DY, all_masks_Jpsi))
                    data = np.concatenate((data_DY_dpToverGenpT, data_Jpsi_dpToverGenpT))
                    stats_ls_pT = get_stats_1Dhist(np.concatenate((data_DY_pT, data_Jpsi_pT)))
                    stats_ls_qd0 = get_stats_1Dhist(np.concatenate((data_DY_qd0, data_Jpsi_qd0)))
                else: 
                    raise KeyError(f"sample_mod={sample_mod} not recognized.")

                # Store dpT/pT data for RooFit's unbinned fit.
                data_dict[eta_key][pT_key][sample_mod][qd0_key]["dpToverpT_vals"] = data
                data_dict[eta_key][pT_key][sample_mod][qd0_key]["stats_ls_pT"] = stats_ls_pT
                data_dict[eta_key][pT_key][sample_mod][qd0_key]["stats_ls_qd0"] = stats_ls_qd0

                num_passed = len(data)
                sum_mask = np.sum(mask)
                assert sum_mask == num_passed
            # End qd0_ls per sample loop.
        # End sample loop.
    # End pT loop.

with open(fullpath_kinem_arrays_pkl,'wb') as outfile:
    pickle.dump(data_dict, outfile, pickle.HIGHEST_PROTOCOL)
print(f"[INFO] data_dict written to pickle file:\n{fullpath_kinem_arrays_pkl}\n")
    
total_muons_original = vdf_concat_MC_2017_DY.count() + vdf_concat_MC_2017_Jpsi.count()
total_muons_found = total_num_muons_DY + total_num_muons_Jpsi
perdif = perc_diff(total_muons_found, total_muons_original)
print(f"Total muons in files: {total_muons_original}")
print(f"Total muons found across all bins: {total_muons_found} (perc. diff. = {perdif:.2f}%)")