"""
# PURPOSE: 
#  This script makes 1 PDF of kinematic distributions in various eta, pT bins. 
#  Eta, pT, and qd0 cuts are applied. 
#  Performs unbinned iterative gaus fits on each distribution. Does not do binned fits yet.
#  Finally a 2D plot of eta vs. pT is made with the final sigma fit in each square.
# NOTES: 
#  Each page of the PDF shows N kinematic distributions, all of which are in the same eta bin. 
#  Each distribution on a page has a different pT bin. 
#  Only runs over 1 year at a time. 
#  Runs on samples made by the Sample class. 
"""
import time
import os
import sys
import math
import pickle
import matplotlib

import ROOT as r
import numpy as np
import matplotlib.pyplot as plt

from Utils_vaex.vaex_fns import vaex_apply_masks, prepare_vaex_df
from Samples.sample_info import Sample
from d0_Studies.KinBin_Info.kinematic_bins import (equal_entry_bin_edges_eta_mod1_wholenum,
                                                   bin_edges_pT_sevenfifths_to1000GeV_wholenum)
from d0_Utils.d0_dicts import label_LaTeX_dict
from d0_Utils.d0_fns import make_binning_array

from Utils_Python.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from Utils_Python.Utils_Plotting import hist_y_label, make_1D_dist, ncolsrows_from_nplots
from Utils_Python.Utils_StatsAndFits import iterative_fit_gaus_unbinned

from matplotlib.backends.backend_pdf import PdfPages

# Hopefully avoids: Tcl_AsyncDelete: async handler deleted by the wrong thread
matplotlib.use('Agg')
#---------------------------#
#----- User Parameters -----#
#---------------------------#
FIXME: outdir = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/dpToverpT/"
filename_base = "gausiterfitsigmas_final"
sample_ls = [
    Sample("MC", "2018", "Jpsi", "hdf5"),
    Sample("MC", "2018", "DY", "hdf5"),
]   
verbose = False 
overwrite = False
iter_gaus = (True, 5)

# Binning.
eta_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40] #equal_entry_bin_edges_eta_mod1_wholenum
pT_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0, 1000.0] #bin_edges_pT_sevenfifths_to1000GeV_wholenum
qd0_range = [-0.020, 0.020]
# qd0_bin_info = [-0.020, 0.020, 0.0002]  # min, max, bin_width.
x_range = [-0.25, 0.25]  # For viewing plots. 
x_bin_info = [-0.25, 0.25, 0.001]

kinem = "delta_pToverRecpT"
title_2D_plot = "Unbinned Iter. Fit #sigma_{Gaus}(#Deltap_{T}/p_{T} dist.) (%)"
color_lim = [0, 12]  # Can also be: None
#-------------------------------------#
#----- Script-Specific Functions -----#
#-------------------------------------#
def check_sample_ls(sample_ls):
    yr_ls = []
    type_ls = []
    for smpl in sample_ls:
        assert smpl.name in ["DY", "Jpsi"]
        yr_ls.append(smpl.year)
        type_ls.append(smpl.data_type)
    assert len(set(yr_ls)) == 1
    assert len(set(type_ls)) == 1

def make_prefix(sample_ls):
    name_str = ""
    for smpl in sample_ls:
        name_str += smpl.name
    return f"{sample_ls[0].data_type}{sample_ls[0].year}{name_str}"

def make_suffix(kinem, num_fits, eta_ls, pT_ls):
    suffix = (
        f"{kinem}_{num_fits}unbinnedfits"
        f"_{min(eta_ls)}eta{max(eta_ls)}"
        f"_{min(pT_ls):.1f}pT{max(pT_ls):.1f}GeV"
        )
    suffix = make_str_title_friendly(suffix)
    return suffix

def make_full_filename(kinem, sample_ls, num_fits, eta_ls, pT_ls, filename_base):
    prefix = make_prefix(sample_ls)
    suffix = make_suffix(kinem, num_fits, eta_ls, pT_ls)
    return f"{prefix}_{filename_base}_{suffix}.pdf"

def prep_area(outdir, outfile_name, overwrite):
    outfile_path = os.path.join(outdir, outfile_name)
    check_overwrite(outfile_path, overwrite=overwrite)
    makeDirs(outdir)
    return outfile_path

def prep_plots(sample_ls, eta_ls, pT_ls):
    n_pages = len(eta_ls)-1
    n_plots = len(pT_ls) - 1
    rows, cols = ncolsrows_from_nplots(n_plots)
    print(f"Making a {n_pages}-page PDF.")
    print(f"Making {n_plots} plots per page ({rows} x {cols} grid).\n")
    print(f"|eta| regions:\n{np.round(eta_ls, decimals=2)}")
    print(f"pT regions:\n{np.round(pT_ls, decimals=2)}\n")
    for smpl in sample_ls:
        print(f"   ...Running over {smpl.data_type} {smpl.year} {smpl.name} sample...")
    return n_pages, n_plots, rows, cols

def make_cut_str(sample_ls, eta_range, pT_range, qd0_range):
    """Return a string of cuts to show on plot."""
    eta_min, eta_max = eta_range[0], eta_range[1]
    pT_min, pT_max = pT_range[0], pT_range[1]
    qd0_min, qd0_max = qd0_range[0], qd0_range[1]

    cuts  = r"$%.2f < \left| \eta^{\mathrm{REC}} \right| < %.2f$," % (eta_min, eta_max) + "\n"
    cuts += r"$%.1f < p_{T} < %.1f$ GeV," % (pT_min, pT_max) + "\n"
    cuts += r"$%.4f < q(\mu)*d_{0}^{\mathrm{%s}} < %.4f$, " % (qd0_min, "BS", qd0_max) + "\n"

    for smpl in sample_ls:
        # Sample-specific items.
        LaTeX_str = smpl.LaTeX_name
        m_inv_cut_min = smpl.inv_mass_cut_lim[0]
        m_inv_cut_max = smpl.inv_mass_cut_lim[1]
        dR_max = smpl.dR_cut
        cuts += r"$%.1f < m_{%s} < %.1f$ GeV,  " % (m_inv_cut_min, LaTeX_str, m_inv_cut_max)
        cuts += r"$\Delta R < %.3f$" % (dR_max) + "\n"
    return cuts

def make_sample_masks(sample_ls, eta_range, pT_range, qd0_range):
    """
    Create masks for eta, pT, qd0 cuts. 
    Also create masks for m_mumu and dR cuts based on sample name.
    """
    for smpl in sample_ls:
        massZ_minmax = smpl.inv_mass_cut_lim
        dR_cut = smpl.dR_cut
        smpl.mask = vaex_apply_masks(smpl.vdf_prepped, eta_range, pT_range, qd0_range, massZ_minmax, dR_cut)

def apply_masks(sample_ls, kinem):
    """
    Make a giant array of masked data from all samples to be fed into a histogram.
    
    kinem must be a column name in a prepped vaex DataFrame.
    """
    data = []
    for smpl in sample_ls:
        vals = smpl.vdf_prepped[kinem].values  # Get col as numpy.ndarray.
        vals_after_cut = vals[smpl.mask.values]  # Apply mask. 
        data.extend(vals_after_cut)
    return np.asarray(data)

def make_ax_labels(kinem, bin_width):
    """
    Make str labels for x-axis and y-axis.
    kinem must be a key in label_LaTeX_dict.
    """
    x_label = label_LaTeX_dict[kinem]["independent_label"]
    x_units = label_LaTeX_dict[kinem]["units"]
    y_label = hist_y_label(bin_width, x_units)
    if len(x_units) > 0:
        x_label += " [{}]".format(x_units)
    return x_label, y_label

def make_2D_hist(kinem, title_2D_plot, eta_ls, pT_ls):
    """
    Make 2D hist which will hold the best-fit sigma values
    from iterated gaus fits of the kinematic distributions. 
    
    The x-bin edges are the pT bin edges.
    The y-bin edges are the eta bin edges.
    """
    hist_2D = r.TH2F("h_2d", title_2D_plot, 
              len(pT_ls)-1, np.array(pT_ls, dtype=float), 
              len(eta_ls)-1, np.array(eta_ls, dtype=float) )

    r.gStyle.SetOptStat(0)

    hist_2D.SetXTitle("p_{T} [GeV]")
    hist_2D.SetYTitle("#left|#eta#right|")

    hist_2D.GetXaxis().SetTitleOffset(1.3)
    hist_2D.SetContour(100)
    if "qd0" in kinem:
        r.gStyle.SetPaintTextFormat("6.0f")  # Number format. Don't measure sub-micron.
    elif "delta_pT" in kinem:
        r.gStyle.SetPaintTextFormat("6.2f")  # Number format.
    return hist_2D

def make_pdf(kinem, hist, outfile_path, preview_plots=False):
    r.gROOT.SetBatch(not preview_plots)
    outf = outfile_path.replace(".pdf", f"_{kinem}bestfitsigma.pdf")
    c1 = r.TCanvas()
    c1.SetLogx(True)
    hist.Draw("colz text")
    c1.SaveAs(outf)

def shift_sigma_val(kinem, bestfit_sigma, bestfit_sigma_err):
    """Scale the sigma and sigma_err up or down depending on the typical value of kinem."""
    if "qd0" in kinem:
        # Factor of 10000 to convert cm to um. 
        bestfit_sigma *= 10000.
        bestfit_sigma_err *= 10000.
    elif "delta_pT" in kinem:
        bestfit_sigma *= 100.  # Convert sigma to percent.
        bestfit_sigma_err *= 100.
    return (bestfit_sigma, bestfit_sigma_err)

def save_stats_to_main_dict(kinem, all_stats_dict, eta_range, pT_range, bestfit_sigma, bestfit_sigma_err):
    key = f"{eta_range[0]}eta{eta_range[1]}_{pT_range[0]}pT{pT_range[1]}_sigmawitherr_{kinem}"
    key = make_str_title_friendly(key)
    all_stats_dict[key] = [bestfit_sigma, bestfit_sigma_err]

def save_dict_to_pkl(all_stats_dict, outpath):
    with open(outpath, "wb") as outf:
        pickle.dump(all_stats_dict, outf, protocol=2)

if __name__ == "__main__":
    # Prep work.
    check_sample_ls(sample_ls)
    outfile_name = make_full_filename(kinem, sample_ls, iter_gaus[1], eta_ls, pT_ls, filename_base)
    outfile_path = prep_area(outdir, outfile_name, overwrite)

    # Plot info.
    plt.style.use("/ufrc/avery/rosedj1/HiggsMassMeasurement/Utils_Python/Plot_Styles_Matplotlib/grid_multiple_plots.mplstyle")
    n_pages, n_plots, rows, cols = prep_plots(sample_ls, eta_ls, pT_ls)
    x_bins, bin_width = make_binning_array(x_bin_info)
    
    # Analysis.
    all_stats_dict = {
        "binedges_eta" : eta_ls, 
        "binedges_pT" : pT_ls
        }
    h_2d_sigma = make_2D_hist(kinem, title_2D_plot, eta_ls, pT_ls)
    # Make only 1 PDF. Each page is a different eta cut.
    with PdfPages(outfile_path) as pdf:
        total_entries = 0
        # Loop over eta regions. 
        for page, (eta_min, eta_max) in enumerate(zip(eta_ls[:-1], eta_ls[1:]), 1):
            eta_range = [eta_min, eta_max]
            
            # Within this eta region, scan the pT regions. 
            t_start = time.perf_counter()
            f = plt.figure()
            for count, (pT_min,pT_max) in enumerate(zip(pT_ls[:-1], pT_ls[1:])):
                pT_range = [pT_min, pT_max]
                ax = plt.subplot(rows,cols,count+1)
                x_label, y_label = make_ax_labels(kinem + "1", bin_width)  # "qd0BS1"

                make_sample_masks(sample_ls, eta_range, pT_range, qd0_range)
                data = apply_masks(sample_ls, kinem)
                ax, bin_vals, bin_edges, stats = make_1D_dist(ax=ax, data=data,
                                                              x_limits=x_range, x_bin_edges=x_bins, 
                                                              x_label=x_label, y_label=y_label,
                                                              title="",
                                                              y_max=-1,
                                                              log_scale=False, color=None, leg_loc=None, 
                                                              framealpha=0.9, display="sci")
                
                cut_str = make_cut_str(sample_ls, eta_range, pT_range, qd0_range)
                ax.text(0.025, 0.78, cut_str, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes,
                            bbox=dict(boxstyle='square', facecolor='white', alpha=0.9))
                
                t_end = time.perf_counter()
                total_entries += stats[0]

                # Optional iterative Gaussian fits.
                if (iter_gaus[0]):
                    # Do iterative fitting procedure.
                    num_iters = iter_gaus[1]
                    # Produce a dictionary of stats for these fits.

                    # if (do_unbinned_fit):
                    # Returns a dictionary with different keys from the binned fit. 
                    fit_stats_dict, ax = iterative_fit_gaus_unbinned(num_iters, data,
                                                                    bin_edges=bin_edges, 
                                                                    bin_vals=bin_vals,
                                                                    num_sigmas=2,
                                                                    ax=ax, draw_on_axes=True, framealpha=0.9, verbose=verbose)
                    #--- It should be easy to implement the else block below.
                    #--- For now I just need to use the unbinned fit. 
                    # else: 
                    #     # Do binned fit.
                    #     # The Gaus coeff is usually around the max height of Gaus curve.
                    #     best_guess_coeff = np.max(bin_vals)  
                    #     best_guess_mean = stats_binned[1]
                    #     best_guess_stdev = stats_binned[3]

                    #     fit_stats_dict, ax = iterative_fit_gaus(num_iters, 
                    #                                         bin_edges, bin_vals, 
                    #                                     param_guess=[best_guess_coeff, 
                    #                                                 best_guess_mean, 
                    #                                                 best_guess_stdev],
                    #                                     param_bounds=([0,-1000,-100], [999999999,1000,100]),
                    #                                     ax=ax, draw_on_axes=True, verbose=False, skip_bad_fit=skip_bad_fit)
                else:
                    fit_stats_dict = None

                bestfit_sigma = fit_stats_dict["stdev_ls"][-1]
                bestfit_sigma_err = fit_stats_dict["stdev_err_ls"][-1]
                pT_midpt = sum(pT_range) / 2.
                eta_midpt = sum(eta_range) / 2.

                shifted_sigma, shifted_sigma_err = shift_sigma_val(kinem, bestfit_sigma, bestfit_sigma_err)
                h_2d_sigma.Fill(pT_midpt, eta_midpt, shifted_sigma)  # Fill(x, y, weight)

                # Save stats.
                save_stats_to_main_dict(kinem, all_stats_dict, eta_range, pT_range, bestfit_sigma, bestfit_sigma_err)
            # End pT loop.
            
            plt.tight_layout()
            pdf.savefig()
            plt.close("all")
            print(f"Page {page}/{n_pages} made. Time taken: {t_end - t_start:.2f} s")
            # del f

        # End eta loop.
    # Close pdf.
    print(f"PDF made at:\n  {outfile_path}")
    print(f"Making PDF of 2D hist...")
    if color_lim is not None:
        h_2d_sigma.GetZaxis().SetRangeUser(*color_lim)
    make_pdf(kinem, h_2d_sigma, outfile_path)

    total_muons_original = 0
    for smpl in sample_ls:
        total_muons_original += smpl.vdf_prepped.count()
    print(f"Total muons in all samples: {total_muons_original}")
    print(f"Total muons found after cuts: {total_entries}")

    outpkl_path = outfile_path.replace(".pdf", ".pkl")
    print(f"Saving stats dict to pickle:\n  {outpkl_path}")
    save_dict_to_pkl(all_stats_dict, outpkl_path)
# End main.