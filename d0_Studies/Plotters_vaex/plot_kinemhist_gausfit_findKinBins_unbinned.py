"""
# WARNING: This script has most likely been superceded by "findKinBins_doFits_makePDFs.py"
#   Do a comparison before deleting this script. -- Jake, 2020-06-24
# PURPOSE: 
#   Make PDFs of any kinematic distribution (e.g., delta_pT/pT), 
#   in specified eta, pT, and q*d0 bins.
#   Also performs unbinned Gaussian fits of distributions, can show the fit stats
#   on the plots, and saves the fit info in a ".stat" file.
#   Optional: make a pickle file of the list of KinBin3D objects. 
# NOTES:   
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
# SYNTAX:  python script.py
# AUTHOR:  Jake Rosenzweig
# UPDATED: 2020-06-24
"""
import os
import sys
import math
import time
import pickle

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Local imports.
from Utils_vaex.vaex_fns import vaex_apply_masks, prepare_vaex_df
from Samples.sample_info import Sample
from d0_Studies.KinBin_Info.kinematic_bins (equal_entry_bin_edges_eta_mod1,
                                        equal_entry_bin_edges_eta_mod2,
                                        equal_entry_bin_edges_eta_mod1_wholenum,
                                        equal_entry_bin_edges_pT_sevenfifths_to1000GeV,
                                        bin_edges_pT_sevenfifths_to1000GeV_wholenum,
                                        binedges_qd0_equalentry_smallest_qd0RMS,
                                        binedges_qd0_equalentry_smallest_qd0RMS_clipped,
                                        binedges_qd0_tracker_res,
                                        binedges_qd0_tracker_res_mod
                                        )
from d0_Utils.d0_dicts import label_LaTeX_dict
from d0_Utils.d0_fns import make_binning_array
from d0_Utils.d0_cls import KinBin3D

from Utils_Python.printing import announce
from Utils_Python.Utils_Physics import perc_diff
from Utils_Python.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from Utils_Python.Utils_Plotting import hist_y_label, make_1D_dist, ncolsrows_from_nplots, get_stats_1Dhist
from Utils_Python.Utils_StatsAndFits import iterative_fit_gaus, iterative_fit_gaus_unbinned

#---------------------------#
#----- User Parameters -----#
#---------------------------#
# Samples with muons treated "individually":
vdf_concat_MC_2018_Jpsi = Sample("MC", "2018", "Jpsi", "hdf5").vdf_prepped
vdf_concat_MC_2018_DY = Sample("MC", "2018", "DY", "hdf5").vdf_prepped
# vdf_concat_MC_2018_DY = #prepare_vaex_df(vdf_MC_2017_DY)
# vdf_concat_MC_2018_Jpsi = #prepare_vaex_df(vdf_MC_2017_Jpsi)

# Dictionary which contains equal-entry q*d0 bin edges.
# !!!!! FIXME: Substitute this with sys.argv[1] or the other one. 

# RUN THIS ONE 
# inpath_3Dbins_pickle_dict = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/kinbin3D_pkls/fullscan_06reg__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"
# inpath_3Dbins_pickle_dict = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/find_best_binning__eta_pT_qd0/Test02_fullscan_06reg_dpToverGenpTsqred__1p4_eta_1p6__30p0_pT_50p0_GeV.pkl"
# inpath_3Dbins_pickle_dict = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/find_best_binning__eta_pT_qd0/fullscan_12reg_3000perreg__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"
# inpath_3Dbins_pickle_dict = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Pickles/fullscan_12reg_3000perreg__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"
# inpath_3Dbins_pickle_dict = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Pickles/Test02_fullscan_06reg_dpToverGenpTsqred__1p4_eta_1p6__30p0_pT_50p0_GeV.pkl"

#--- Delete below:
inpath_3Dbins_pickle_dict = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/20200624_MC2018JpsiDY_test01_6regwith3000perreg__0p0_eta_0p4__5p0_pT_1000p0_GeV.pkl"
#--- Delete above:

# 20200604: run below
# inpath_3Dbins_pickle_dict = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/20200624_fullstats_MC2018JpsiDY_6regwith3000perreg__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"

#-----Below is a big boy. -----#
# inpath_3Dbins_pickle_dict = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/find_best_binning__eta_pT_qd0/equalentry_qd0_binedges__20_regions_max__atleast1000entriesperregion__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"
# outdir_plots = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/hists_dpToverpT/MC/2017/"
outdir_plots = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/hists_dpToverpT/"
# outdir_kinbin_pkl = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/kinbin3D_pkls/"
outdir_kinbin_pkl = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/kinbins/"

# filename_base = "realtest01_MC2017DYandJpsi_fullscan__autodetect_qd0_regions"
# filename_base = "fullpath_fullscan_12reg"
filename_base = "20200624_testscan_06reg_test01"
write_to_pickle = True
overwrite = False
verbose = True
skip_bad_fit = True
do_unbinned_fit = True  # If false, will default to binned fit.
iter_gaus = (True, 5)

# Binning.
# eta_ls = equal_entry_bin_edges_eta_mod1_wholenum
# pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum
eta_ls = [0.0, 0.4]
pT_ls = [20.0, 27.0, 38.0, 50.0]

# Kinematic to be plotted on all histograms. 
# Should not contain 1 or 2. 
# Acceptable values found in prepare_vaex_df().
kinem = "delta_pToverGenpT"  
x_bin_info = [-0.3, 0.3, 0.003]
x_zoom_range = [-0.35, 0.35]

# Cuts to make.
dR_max_DY = 0.002
dR_max_Jpsi = 0.005
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
makeDirs(outdir_plots)

suffix   = "_{}_eta_{}".format(min(eta_ls), max(eta_ls))
suffix  += "_{:.1f}_pT_{:.1f}_GeV".format(min(pT_ls), max(pT_ls))
if (do_unbinned_fit): 
    suffix += "_unbinned"
suffix = make_str_title_friendly(suffix)

stats_fullpath = os.path.join(outdir_plots, filename_base + suffix + ".stat")
fullpath_kinbin_ls_pkl = os.path.join(outdir_kinbin_pkl, filename_base + suffix + ".pkl") 

check_overwrite(stats_fullpath, overwrite=overwrite) 
check_overwrite(fullpath_kinbin_ls_pkl, overwrite=overwrite) 
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

with open(inpath_3Dbins_pickle_dict, "rb") as f:
    equalentry_binedge_dict = pickle.load(f)

# plt.style.use('grid_multiple_plots')
plt.style.use("/home/rosedj1/.config/matplotlib/mpl_configdir/stylelib/grid_multiple_plots_python2.mplstyle")

total_entries = 0
kinbin3D_ls = []

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

# Unpack kinematic data.
eta_arr_DY = vdf_concat_MC_2018_DY.evaluate("eta")
pT_arr_DY = vdf_concat_MC_2018_DY.evaluate("pT")
qd0_arr_DY = vdf_concat_MC_2018_DY.evaluate("qd0BS")
massZ_arr_DY = vdf_concat_MC_2018_DY.evaluate("massZ")
dR_arr_DY = vdf_concat_MC_2018_DY.evaluate("delta_R")

eta_arr_Jpsi = vdf_concat_MC_2018_Jpsi.evaluate("eta")
pT_arr_Jpsi = vdf_concat_MC_2018_Jpsi.evaluate("pT")
qd0_arr_Jpsi = vdf_concat_MC_2018_Jpsi.evaluate("qd0BS")
massZ_arr_Jpsi = vdf_concat_MC_2018_Jpsi.evaluate("massZ")
dR_arr_Jpsi = vdf_concat_MC_2018_Jpsi.evaluate("delta_R")

# Kinematic to be plotted in each KinBin3D.
kinem_arr_DY = vdf_concat_MC_2018_DY.evaluate(kinem)
kinem_arr_Jpsi = vdf_concat_MC_2018_Jpsi.evaluate(kinem)
#----------------#
#----- Main -----#
#----------------#

# Loop over eta regions.
for k in range(len(eta_ls)-1):
    eta_min = eta_ls[k]
    eta_max = eta_ls[k+1]
    eta_range = [eta_min, eta_max]
    print("eta loop k,",k)

    # For each eta range, make a pdf. 
    extra   = "__{}_eta_{}".format(eta_min, eta_max)
    extra  += "__{:.1f}_pT_{:.1f}_GeV".format(min(pT_ls), max(pT_ls))
    extra = make_str_title_friendly(extra)
    extra += ".pdf"

    outpath_pdf = os.path.join(outdir_plots, filename_base + extra)
    check_overwrite(outpath_pdf, overwrite=overwrite)    

    print("Making PDF {}/{}...".format(k+1, n_pdfs))
    t_start_page = time.perf_counter()
    with PdfPages(outpath_pdf) as pdf:
        # Within 1 eta region, scan the pT regions. 
        for j in range(len(pT_ls)-1):
            pT_min = pT_ls[j]
            pT_max = pT_ls[j+1]
            pT_range = [pT_min, pT_max]
            print("pT loop j,",j)

            # Within 1 pT region, scan the q*d0 regions and make 1 page of PDF. 
            t_start = time.perf_counter()
            f = plt.figure()

            eta_key = "eta_bin_left_edge={}".format(eta_min)
            pT_key = "pT_bin_left_edge={}".format(pT_min)
            qd0_ls = equalentry_binedge_dict[eta_key][pT_key]
            print("eta_key,",eta_key)
            print("pT_key,",pT_key)
            print("qd0_ls,",qd0_ls)

            rows, cols = get_grid_info(qd0_ls)
            for count in range(len(qd0_ls)-1):
                ax = plt.subplot(rows,cols,count+1)
                qd0_min = qd0_ls[count]
                qd0_max = qd0_ls[count+1]
                qd0_range = [qd0_min, qd0_max]
                print("qd0 loop", count)
                print("qd0_min {}, qd0_max {}, qd0_range {}".format(qd0_min, qd0_max, qd0_range))
                
                x_bins, bin_width = make_binning_array(x_bin_info)

                x_label = label_LaTeX_dict[kinem]["independent_label"]
                x_units = label_LaTeX_dict[kinem]["units"]
                y_label = hist_y_label(bin_width, x_units)
                if len(x_units) > 0:
                    x_label += " [{}]".format(x_units)

                cuts  = r"$%.2f < \left| \eta^{\mathrm{REC}} \right| < %.2f$," % (eta_min, eta_max) + "\n"
                cuts += r"$%.1f <$ %s $< %.1f$ GeV," % (pT_min, p_str_latex, pT_max) + "\n"
                cuts += r"$%.4f < q(\mu)*d_{0}^{\mathrm{%s}} < %.4f$, " % (qd0_min, "BS", qd0_max) + "\n"
                cuts += r"$%.1f < m_{J/\psi} < %.1f$ GeV,  " % (massZ_min_Jpsi, massZ_max_Jpsi)
                cuts += r"$\Delta R < %.3f$" % (dR_max_Jpsi) + "\n"
                cuts += r"$%.1f < m_{Z} < %.1f$ GeV,  " % (massZ_min_DY, massZ_max_DY)
                cuts += r"$\Delta R < %.3f$" % (dR_max_DY)
                print("cuts str,",cuts)
                
                # At this point, vdf_concat_MC_2017_X definitely has full entries.
                # Previous cuts not applied here.
                print("\napply masks with these cuts:")
                print("eta_range",eta_range)
                print("pT_range",pT_range)
                print("qd0_range",qd0_range)
                print("For DY: massZ_minmax_DY",massZ_minmax_DY)
                print("dR_max_DY",dR_max_DY)
                print("For Jpsi: massZ_minmax_Jpsi",massZ_minmax_Jpsi)
                print("dR_max_Jpsi",dR_max_Jpsi)

                # Prepare masks.
                mask_eta_DY = (eta_min < np.abs(eta_arr_DY)) & (np.abs(eta_arr_DY) < eta_max)
                mask_pT_DY = (pT_min < pT_arr_DY) & (pT_arr_DY < pT_max)
                mask_qd0_DY = (qd0_min < qd0_arr_DY) & (qd0_arr_DY < qd0_max)
                mask_massZ_DY = (massZ_min_DY < massZ_arr_DY) & (massZ_arr_DY < massZ_max_DY)
                mask_dR_DY = (dR_arr_DY < dR_max_DY)

                mask_eta_Jpsi = (eta_min < np.abs(eta_arr_Jpsi)) & (np.abs(eta_arr_Jpsi) < eta_max)
                mask_pT_Jpsi = (pT_min < pT_arr_Jpsi) & (pT_arr_Jpsi < pT_max)
                mask_qd0_Jpsi = (qd0_min < qd0_arr_Jpsi) & (qd0_arr_Jpsi < qd0_max)
                mask_massZ_Jpsi = (massZ_min_Jpsi < massZ_arr_Jpsi) & (massZ_arr_Jpsi < massZ_max_Jpsi)
                mask_dR_Jpsi = (dR_arr_Jpsi < dR_max_Jpsi)

                all_masks_DY = mask_eta_DY & mask_pT_DY & mask_qd0_DY & mask_massZ_DY & mask_dR_DY
                all_masks_Jpsi = mask_eta_Jpsi & mask_pT_Jpsi & mask_qd0_Jpsi & mask_massZ_Jpsi & mask_dR_Jpsi

                # The data to be histogrammed.
                selected_muons_DY = kinem_arr_DY[all_masks_DY]
                selected_muons_Jpsi = kinem_arr_Jpsi[all_masks_Jpsi]

                # Get data and stats of q*d0 within this 3Dcube.
                qd0_vals_this3Dcube_DY = qd0_arr_DY[all_masks_DY]
                qd0_vals_this3Dcube_Jpsi = qd0_arr_Jpsi[all_masks_Jpsi]
                qd0_vals_this3Dcube = np.concatenate( (qd0_vals_this3Dcube_DY, qd0_vals_this3Dcube_Jpsi) )
                qd0_this3Dcube_stats_ls = get_stats_1Dhist(qd0_vals_this3Dcube)
                print("qd0_this3Dcube_stats_ls:",qd0_this3Dcube_stats_ls)

                # Get data and stats of pT within this 3Dcube.
                pT_vals_this3Dcube_DY = pT_arr_DY[all_masks_DY]
                pT_vals_this3Dcube_Jpsi = pT_arr_Jpsi[all_masks_Jpsi]
                pT_vals_this3Dcube = np.concatenate( (pT_vals_this3Dcube_DY, pT_vals_this3Dcube_Jpsi) )
                pT_this3Dcube_stats_ls = get_stats_1Dhist(pT_vals_this3Dcube)
                print("pT_this3Dcube_stats_ls:",pT_this3Dcube_stats_ls)

                num_passed_DY = len(selected_muons_DY)
                num_passed_Jpsi = len(selected_muons_Jpsi)
                sum_mask_DY = np.sum(all_masks_DY)
                sum_mask_Jpsi = np.sum(all_masks_Jpsi)
                assert sum_mask_DY == num_passed_DY
                assert sum_mask_Jpsi == num_passed_Jpsi
                assert len(qd0_vals_this3Dcube_DY) == num_passed_DY
                assert len(qd0_vals_this3Dcube_Jpsi) == num_passed_Jpsi

                print("np.sum(all_masks_DY),", sum_mask_DY)
                print("len(selected_muons_DY),", num_passed_DY)
                print("np.sum(all_masks_Jpsi),", sum_mask_Jpsi)
                print("len(selected_muons_Jpsi),", num_passed_Jpsi)

                n_entries_this3Dcube = num_passed_DY + num_passed_Jpsi
                print("According to n_entries_this3Dcube: {}".format(n_entries_this3Dcube))
                total_entries += n_entries_this3Dcube
                print("total_entries cumulative,",total_entries,"\n")

                # selected_DY_arr   = vdf_concat_MC_2018_DY.evaluate(  kinem,selection=all_masks_DY)
                # selected_Jpsi_arr = vdf_concat_MC_2018_Jpsi.evaluate(kinem,selection=all_masks_Jpsi)
                
                # data = np.concatenate((selected_DY_arr, selected_Jpsi_arr))
                data = np.concatenate((selected_muons_DY, selected_muons_Jpsi))

                announce("--- Doing unbinned fit ---")
                # Plot the kinem hist.
                ax, bin_vals, bin_edges, stats = make_1D_dist(
                                                    ax=ax, 
                                                    # data definitely gets cut as expected.
                                                    data=data,
                                                    x_limits=x_zoom_range,
                                                    x_bins=x_bins, 
                                                    x_label=x_label, 
                                                    y_label=y_label,
                                                    title="",
                                                    y_max=-1,
                                                    log_scale=False)
                print("stats ls: {}".format(stats))
                ax.text(0.025, 0.78, cuts, horizontalalignment='left', verticalalignment='center', 
                transform=ax.transAxes, bbox=dict(boxstyle='square', facecolor='white', alpha=0.9))

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
                        # Use plotted kinem as the key for this dict of stats. 
                        # The Gaus coeff is usually around the max height of Gaus curve.
                        best_guess_coeff = np.max(bin_vals)  
                        best_guess_mean = stats[1]
                        best_guess_stdev = stats[3]

                        fit_stats_dict, ax = iterative_fit_gaus(num_iters, 
                                                                bin_edges, bin_vals, 
                                                            param_guess=[best_guess_coeff, 
                                                                        best_guess_mean, 
                                                                        best_guess_stdev],
                                                            param_bounds=([0,-1000,-100], [999999999,1000,100]),
                                                            ax=ax, draw_on_axes=True, verbose=False, skip_bad_fit=skip_bad_fit)
                    
                # Save all them tasty fit statistics. 
                kinbin3D_ls.append( KinBin3D(
                                        eta_range=eta_range,
                                        pT_range=pT_range, 
                                        qd0_range=qd0_range,
                                        n_entries=n_entries_this3Dcube,
                                        kinem=kinem,
                                        fit_stats_dict=fit_stats_dict,
                                        pT_stats_ls=pT_this3Dcube_stats_ls,
                                        qd0_stats_ls=qd0_this3Dcube_stats_ls,
                                    )
                                  )

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

            msg = "Page {}/{} made. Time taken: {:.2f} s".format(j+1, n_pages, t_end - t_start)
            announce(msg)
            print()

        # End pT loop. Make next page.
    t_end_page = time.perf_counter()
    msg_made_pdf = "PDF {}/{} made".format(k+1, n_pdfs)
    announce(msg_made_pdf, pad_char="@", n_center_pad_chars=5)
    print("location:\n  {}".format(outpath_pdf))
    print("(took {:.2f} s)\n".format(t_end_page - t_start_page))
    # Save this 1 eta reg, 1 pT reg, and all q*d0 plots on one page. 

# Go to next pT reg and next page.
print("All PDFs created.")

if (write_to_pickle):
    with open(fullpath_kinbin_ls_pkl,'wb') as output:
        pickle.dump(kinbin3D_ls, output, protocol=2)
    print("[INFO] KinBin3D list written to pickle file:\n{}\n".format(fullpath_kinbin_ls_pkl))
    
total_muons_original = vdf_concat_MC_2018_DY.count() + vdf_concat_MC_2018_Jpsi.count()
print("Total muons: {}".format(total_muons_original))
perdif = perc_diff(total_entries, total_muons_original)
print("Total muons found across all bins: {} (perc. diff. = {:.2f}%)".format(total_entries, perdif))
